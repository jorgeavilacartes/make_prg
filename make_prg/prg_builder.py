import logging
from typing import List, Tuple

from make_prg.from_msa import MSA
from make_prg.io_utils import load_alignment_file
from make_prg.from_msa.cluster_sequences import kmeans_cluster_seqs_in_interval
from make_prg.seq_utils import (
    ambiguous_bases,
    remove_duplicates,  # TODO: quantify and remove dups
    get_interval_seqs,
    NONMATCH,
)
from make_prg.from_msa.interval_partition import IntervalPartitioner, Interval, IntervalType
import pickle
from pathlib import Path
import os
import shlex
import time
from Bio import SeqIO
from abc import ABC, abstractmethod
import pyabpoa as pa
import uuid
import shutil
import subprocess
import copy
import numpy as np
from Bio.Seq import Seq


class MSAAligner(ABC):
    @classmethod
    @abstractmethod
    def get_updated_alignment(cls, leaf_name: str, previous_alignment: Path, new_sequences: List[str], temp_prefix: Path) -> Path:
        pass

class MSAAlignerAbPOA(MSAAligner):
    aligner = pa.msa_aligner(aln_mode='g',
                             extra_b=-1,  # adaptive banding disabled
                             is_diploid=0)

    @classmethod
    def get_updated_alignment(cls, leaf_name: str, previous_alignment: Path, new_sequences: List[str], temp_prefix: Path) -> Path:
        previous_alignment = shlex.quote(str(previous_alignment))
        start = time.time()

        # TODO: AbPOA do not process multiline fasta!!!! Run sample example and see results
        msa_result = cls.aligner.msa(new_sequences, False, True, incr_fn=previous_alignment)
        new_msa_filepath = temp_prefix / "updated_msa.fa"
        with open(new_msa_filepath, "w") as new_msa_filehandler:
            for index_seq, seq in enumerate(msa_result.msa_seq):
                print(f">Denovo_path_{index_seq}_random_id_{uuid.uuid4()}", file=new_msa_filehandler)
                print(seq, file=new_msa_filehandler)

        stop = time.time()
        runtime = stop-start
        print(f"abPOA update runtime for {leaf_name} in seconds: {runtime:.3f}")

        return new_msa_filepath


class MSAAlignerMAFFT(MSAAligner):
    @classmethod
    def get_updated_alignment(cls, leaf_name: str, previous_alignment: Path, new_sequences: List[str], temp_prefix: Path) -> Path:
        new_sequences_filename = temp_prefix / "new_sequences.fa"
        with open(new_sequences_filename, "w") as new_sequences_handler:
            for index_new_seq, new_seq in enumerate(new_sequences):
                print(f">Denovo_path_{index_new_seq}_random_id_{uuid.uuid4()}", file=new_sequences_handler)
                print(new_seq, file=new_sequences_handler)

        new_msa = temp_prefix / "updated_msa.fa"
        previous_alignment = shlex.quote(str(previous_alignment))
        mafft_tmpdir = Path(temp_prefix / f"temp.mafft.{uuid.uuid4()}")
        if mafft_tmpdir.exists():
            shutil.rmtree(mafft_tmpdir)
        mafft_tmpdir.mkdir()
        env = os.environ
        env["TMPDIR"] = str(mafft_tmpdir)

        args = " ".join(
            [
                "mafft",
                "--auto",
                "--quiet",
                "--thread",
                "1",
                "--add",
                str(new_sequences_filename),
                str(previous_alignment),
                ">",
                str(new_msa),
            ]
        )

        start = time.time()
        process = subprocess.Popen(
            args, stderr=subprocess.PIPE, encoding="utf-8", shell=True, env=env,
        )
        exit_code = process.wait()
        shutil.rmtree(mafft_tmpdir)
        if exit_code != 0:
            raise RuntimeError(
                f"Failed to execute MAFFT for {leaf_name} due to the following error:\n"
                f"{process.stderr.read()}"
            )
        stop = time.time()
        runtime = stop-start
        print(f"MAFFT update runtime for {leaf_name} in seconds: {runtime:.3f}")

        return new_msa


class PrgBuilderRecursiveTreeNode(ABC):
    def __init__(self,
                 nesting_level,
                 alignment,
                 parent,
                 prg_builder):
        # set the basic attributes
        self.id = prg_builder.get_next_node_id()
        self.nesting_level = nesting_level
        self.parent = parent
        self.prg_builder = prg_builder
        self.new_sequences = None

        self.alignment = self.remove_gaps(alignment)

        self._set_derived_helper_attributes()

        # generate recursion tree
        self._children = self._get_children()

    def remove_gaps(self, alignment):
        """
        Return a gapless alignment. This code is long and a bit convoluted because it is optimised (it was too slow if
        done in the most intuitive way).
        """
        alignment_as_array = np.array([list(rec) for rec in alignment], str, order="F")
        gapless_sequences = [[] for _ in range(len(alignment))]
        for column_index in range(alignment.get_alignment_length()):
            column_bases = alignment_as_array[:, column_index]
            column_bases_deduplicated = list(set(column_bases))
            just_gaps = column_bases_deduplicated == ["-"]
            if not just_gaps:
                for gapless_sequence, base in zip(gapless_sequences, column_bases):
                    gapless_sequence.append(base)

        gapless_records = []
        for gapless_sequence, previous_record in zip(gapless_sequences, alignment):
            new_record = copy.deepcopy(previous_record)
            new_record.seq = Seq("".join(gapless_sequence))
            gapless_records.append(new_record)

        gapless_alignment = MSA(gapless_records)

        if gapless_alignment.get_alignment_length() == 0:
            print(f"Alignment composed only of gaps: {self.prg_builder.msa_file}")

        return gapless_alignment



    @abstractmethod
    def _set_derived_helper_attributes(self):
        pass

    @abstractmethod
    def _get_children(self):
        pass

    @abstractmethod
    def preorder_traversal_to_build_prg(self, prg_as_list, delim_char):
        pass

    @abstractmethod
    def batch_update(self, temp_prefix):
        pass


class PrgBuilderMultiClusterNode(PrgBuilderRecursiveTreeNode):
    def __init__(self,
                 nesting_level,
                 alignment,
                 parent,
                 prg_builder):
        super().__init__(
            nesting_level,
            alignment,
            parent,
            prg_builder
        )

    def _set_derived_helper_attributes(self):
        pass  # nothing to set here

    def _get_children(self):
        # each child is a PrgBuilderSingleClusterNode for each cluster subalignment
        cluster_subalignments = self._get_subalignments_by_clustering()
        children = []
        for alignment in cluster_subalignments:
            child = PrgBuilderSingleClusterNode(
                nesting_level=self.nesting_level,
                alignment=alignment,
                parent=self,
                prg_builder=self.prg_builder
            )
            children.append(child)
        return children

    ##################################################################################
    # traversal methods
    def preorder_traversal_to_build_prg(self, prg_as_list, delim_char):
        site_num = self.prg_builder.get_next_site_num()
        prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")

        for child_index, child in enumerate(self._children):
            site_num_to_separate_alleles = (site_num + 1) if (child_index < len(self._children) - 1) else site_num
            child.preorder_traversal_to_build_prg(prg_as_list, delim_char)
            prg_as_list.extend(f"{delim_char}{site_num_to_separate_alleles}{delim_char}")
    ##################################################################################


    #####################################################################################################
    #  Clustering methods
    def _get_subalignments_by_clustering(self):
        id_lists = kmeans_cluster_seqs_in_interval(
            self.alignment,
            self.prg_builder.min_match_length,
        )
        list_sub_alignments = [self._get_sub_alignment_by_list_id(id_list) for id_list in id_lists]
        return list_sub_alignments

    def _get_sub_alignment_by_list_id(
        self, id_list: List[str]
    ):
        list_records = [record for record in self.alignment if record.id in id_list]
        sub_alignment = MSA(list_records)
        return sub_alignment
    #####################################################################################################

    def batch_update(self, temp_prefix):
        assert False, "Update was called on a non-leaf node"


class PrgBuilderSingleClusterNode(PrgBuilderRecursiveTreeNode):
    def __init__(self,
                 nesting_level,
                 alignment,
                 parent,
                 prg_builder):
        super().__init__(
            nesting_level,
            alignment,
            parent,
            prg_builder
        )

    def _set_derived_helper_attributes(self):
        self.consensus = self._get_consensus()
        self.length = len(self.consensus)
        (
            self.match_intervals,
            self.non_match_intervals,
            self.all_intervals,
        ) = IntervalPartitioner(
            self.consensus, self.prg_builder.min_match_length, self.alignment
        ).get_intervals()

    def _get_children(self):
        # stop conditions
        single_match_interval = (len(self.all_intervals) == 1) and (self.all_intervals[0] in self.match_intervals)
        max_nesting_level_reached = self.nesting_level == self.prg_builder.max_nesting
        small_variant_site = self.alignment.get_alignment_length() < self.prg_builder.min_match_length
        if single_match_interval or max_nesting_level_reached or small_variant_site:
            return list()

        children = []
        for interval in self.all_intervals:
            sub_alignment = self.alignment[:, interval.start: interval.stop + 1]

            if interval in self.match_intervals:
                # all seqs are not necessarily exactly the same: some can have 'N'
                # thus still process all of them, to get the one with no 'N'.
                seqs = get_interval_seqs(sub_alignment)
                assert len(seqs) == 1, "Got >1 filtered sequences in match interval"
                subclass = PrgBuilderSingleClusterNode
            else:
                subclass = PrgBuilderMultiClusterNode

            child = subclass(
                nesting_level=self.nesting_level + 1,
                alignment=sub_alignment,
                parent=self,
                prg_builder=self.prg_builder
            )
            children.append(child)

        return children


    ##################################################################################
    # properties
    @property
    def prop_in_match_intervals(self):
        length_match_intervals = 0
        for interval in self.match_intervals:
            length_match_intervals += interval.stop - interval.start + 1
        return length_match_intervals / float(self.length)

    @property
    def num_seqs(self):
        return len(self.alignment)
    ##################################################################################


    ##################################################################################
    # traversal methods
    def preorder_traversal_to_build_prg(self, prg_as_list, delim_char):
        is_leaf_node = len(self._children) == 0
        if is_leaf_node:
            self._get_prg(prg_as_list, delim_char)
        else:
            for child in self._children:
                child.preorder_traversal_to_build_prg(prg_as_list, delim_char)
    ##################################################################################

    ##################################################################################
    # helpers
    def _get_consensus(self):
        """ Produces a 'consensus string' from an MSA: at each position of the
        MSA, the string has a base if all aligned sequences agree, and a "*" if not.
        IUPAC ambiguous bases result in non-consensus and are later expanded in the prg.
        N results in consensus at that position unless they are all N."""
        consensus_string = ""
        for i in range(self.alignment.get_alignment_length()):
            column = set([record.seq[i] for record in self.alignment])
            column = column.difference({"N"})
            if (
                len(ambiguous_bases.intersection(column)) > 0
                or len(column) != 1
            ):
                consensus_string += NONMATCH
            else:
                consensus_string += column.pop()
        return consensus_string

    def _get_prg(self, prg_as_list, delim_char):
        for interval in self.all_intervals:
            sub_alignment = self.alignment[:, interval.start : interval.stop + 1]
            seqs = get_interval_seqs(sub_alignment)

            single_seq = len(seqs) == 1
            if single_seq:
                start_index = len(prg_as_list)
                prg_as_list.extend(seqs[0])
                end_index = len(prg_as_list)+1
                self.prg_builder.update_leaves_index(start_index, end_index, node=self)
            else:
                # Add the variant seqs to the prg.
                site_num = self.prg_builder.get_next_site_num()
                prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")
                for seq_index, seq in enumerate(seqs):
                    # TODO: recheck this
                    site_num_for_this_seq = (site_num + 1) if (seq_index < len(seqs) - 1) else site_num
                    start_index = len(prg_as_list)
                    prg_as_list.extend(seq)
                    end_index = len(prg_as_list) + 1
                    self.prg_builder.update_leaves_index(start_index, end_index, node=self)
                    prg_as_list.extend(f"{delim_char}{site_num_for_this_seq}{delim_char}")
    ##################################################################################

    ##################################################################################
    # update methods
    def add_seq_to_batch_update(self, new_sequence: str):
        if self.new_sequences is None:
            self.new_sequences = set()
        self.new_sequences.add(new_sequence)


    def batch_update(self, temp_prefix):
        no_update_to_be_done = self.new_sequences is None or len(self.new_sequences) == 0
        if no_update_to_be_done:
            return
        self._update_leaf(temp_prefix)


    def _update_leaf(self, temp_prefix):
        # update the MSA
        temp_prefix = temp_prefix / self.prg_builder.locus_name / f"node_{self.id}"
        os.makedirs(temp_prefix, exist_ok=True)

        previous_msa_filename = temp_prefix / "previous_msa.fa"
        with open(previous_msa_filename, "w") as previous_msa_handler:
            SeqIO.write(self.alignment, previous_msa_handler, "fasta")

        print(f"Updating MSA for {self.prg_builder.locus_name}, node {self.id}...")
        msa_aligner = MSAAlignerMAFFT
        new_msa = msa_aligner.get_updated_alignment(leaf_name=f"{self.prg_builder.locus_name}, node {self.id}",
                                                   previous_alignment=previous_msa_filename,
                                                   new_sequences=list(self.new_sequences),
                                                   temp_prefix=temp_prefix)

        # update the alignment
        self.alignment = load_alignment_file(str(new_msa), "fasta")

        # update the other fields
        self._set_derived_helper_attributes()

        # reset the new sequences
        self.new_sequences = None

        # regenerate recursion tree
        self._children = self._get_children()
    ##################################################################################


class PrgBuilder(object):
    """
    Prg builder based from a multiple sequence alignment.
    Note min_match_length must be strictly greater than max_nesting + 1.
    """

    def __init__(
        self,
        locus_name,
        msa_file,
        alignment_format,
        max_nesting,
        min_match_length
    ):
        self.locus_name = locus_name
        self.msa_file = msa_file
        self.alignment_format = alignment_format
        self.max_nesting = max_nesting
        self.min_match_length = min_match_length
        self.leaves_index = {}
        self.node_id = 0

        alignment = load_alignment_file(msa_file, alignment_format)
        self._root = PrgBuilderSingleClusterNode(nesting_level=1,
                                                 alignment=alignment,
                                                 parent=None,
                                                 prg_builder=self)

    def build_prg(self):
        self._site_num = 5
        prg_as_list = []
        self._root.preorder_traversal_to_build_prg(prg_as_list, delim_char=" ")
        prg = "".join(prg_as_list)
        return prg

    def get_next_site_num(self):
        previous_site_num = self._site_num
        self._site_num+=2
        return previous_site_num

    def get_next_node_id(self):
        self.node_id += 1
        return self.node_id - 1

    def update_leaves_index(self, start_index: int, end_index: int, node):
        interval = (start_index, end_index)
        self.leaves_index[interval] = node

    def get_node_given_interval(self, interval: Tuple[int, int]):
        # TODO: move this back to assert
        interval_is_indexed = interval in self.leaves_index
        if not interval_is_indexed:
            raise RuntimeError(f"Queried interval {interval} does not exist in leaves index.\n"
                               f"self.locus_name = {self.locus_name}\n"
                               f"self.leaves_index.keys() = {self.leaves_index.keys()}")

        # assert interval in self.leaves_index, \
        #     f"Fatal error: queried interval {interval} does not exist in leaves index"

        return self.leaves_index[interval]

    def serialize(self, filename):
        with open(filename, "wb") as filehandler:
            pickle.dump(self, filehandler)

    @staticmethod
    def deserialize(filename):
        with open(filename, "rb") as filehandler:
            return pickle.load(filehandler)


class PrgBuilderCollection:
    """
    Represent a collection of PrgBuilder and some other info, to be serialised and deserialised
    """
    def __init__(self, locus_name_to_pickle_files, cl_options):
        self.locus_name_to_pickle_files = locus_name_to_pickle_files
        self.cl_options = cl_options

    def serialize(self):
        with open(f"{self.cl_options.output_prefix}.update_DS", "wb") as filehandler:
            pickle.dump(self, filehandler)

    @staticmethod
    def deserialize(filename):
        with open(filename, "rb") as filehandler:
            return pickle.load(filehandler)

    def to_absolute_paths(self):
        output_prefix_parent = Path(self.cl_options.output_prefix).parent.absolute()
        for locus_name, pickle_file in self.locus_name_to_pickle_files.items():
            absolute_filepath = output_prefix_parent / pickle_file
            self.locus_name_to_pickle_files[locus_name] = str(absolute_filepath)