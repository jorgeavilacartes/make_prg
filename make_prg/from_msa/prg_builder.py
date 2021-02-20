import logging
from typing import List

from make_prg.from_msa import MSA
from make_prg.io_utils import load_alignment_file
from make_prg.from_msa.cluster_sequences import kmeans_cluster_seqs_in_interval
from make_prg.seq_utils import (
    ambiguous_bases,
    remove_duplicates,
    get_interval_seqs,
    NONMATCH,
)
from make_prg.from_msa.interval_partition import IntervalPartitioner, Interval, IntervalType
import pickle
from pathlib import Path
import uuid
import shutil
import os
import shlex
import subprocess
import time
from Bio.Seq import Seq
from Bio import SeqIO
from abc import ABC, abstractmethod


# TODO: change to spoa
class MSAAligner:
    @classmethod
    def get_updated_alignment(cls, locus_name: str, previous_alignment: Path, new_sequences: Path, prefix: str) -> Path:
        logging.info(f"Updating MSA for {locus_name}...")

        # make paths shell-safe
        new_msa = shlex.quote(f"{prefix}.updated_msa.fa")
        previous_alignment = shlex.quote(str(previous_alignment))
        new_sequences = shlex.quote(str(new_sequences))
        mafft_tmpdir = Path(f"{prefix}.mafft.{uuid.uuid4()}")
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
                new_sequences,
                previous_alignment,
                ">",
                new_msa,
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
                f"Failed to execute mafft for {locus_name} due to the following error:\n"
                f"{process.stderr.read()}"
            )
        stop = time.time()
        runtime = stop-start
        logging.info(f"Finished updating MSA for {locus_name}")
        logging.info(f"MAFFT update runtime for {locus_name} in seconds: {runtime:.3f}")

        return Path(new_msa)


class PrgBuilderRecursiveTreeNode(ABC):
    def __init__(self,
                 nesting_level,
                 alignment,
                 interval,
                 parent,
                 prg_builder):
        # set the basic attributes
        self.nesting_level = nesting_level
        self.alignment = alignment
        self.interval = interval
        self.parent = parent
        self.prg_builder = prg_builder
        self.new_sequences = None

        self._set_derived_helper_attributes()

        # generate recursion tree
        self._children = self._get_children()

    @abstractmethod
    def _set_derived_helper_attributes(self):
        pass

    @abstractmethod
    def _get_children(self):
        pass

    @abstractmethod
    def preorder_traversal_to_build_prg(self, prg_as_list, delim_char):
        pass


class PrgBuilderMultiClusterNode(PrgBuilderRecursiveTreeNode):
    def __init__(self,
                 nesting_level,
                 alignment,
                 interval,
                 parent,
                 prg_builder):
        super().__init__(
            nesting_level,
            alignment,
            interval,
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
                interval=self.interval,
                parent=self,
                prg_builder=self.prg_builder
            )
            children.append(child)
        return children

    ##################################################################################
    # traversal methods
    # TODO memoize this?
    def preorder_traversal_to_build_prg(self, prg_as_list, delim_char):
        site_num = self.prg_builder.get_next_site_num()
        prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")

        for child_index, child in enumerate(self._children):
            site_num_to_separate_alleles = (site_num + 1) if (child_index < len(self._children) - 1) else site_num
            child.preorder_traversal_to_build_prg(prg_as_list, delim_char)
            prg_as_list.extend(f"{delim_char}{site_num_to_separate_alleles}{delim_char}")
    ##################################################################################

    ##################################################################################
    # update methods
    # TODO: ignore for now
    def add_seqs_to_batch_update(self, new_sequences: List[Seq]):
        if self.new_sequences is None:
            self.new_sequences = []
        self.new_sequences.extend(new_sequences)
        self.new_sequences = list(set(self.new_sequences)) # dedup

    def batch_update(self):
        no_update_to_be_done = self.new_sequences is None or len(self.new_sequences) == 0
        if no_update_to_be_done:
            return

        # update the MSA
        previous_msa_filename = f"{self.prg_builder.prefix}.previous_msa.fa"
        with open(previous_msa_filename, "w") as previous_msa_handler:
            SeqIO.write(self.alignment, previous_msa_handler, "fasta")
        new_sequences_filename = f"{self.prg_builder.prefix}.new_sequences.fa"
        with open(new_sequences_filename, "w") as new_sequences_handler:
            SeqIO.write(self.new_sequences, new_sequences_handler, "fasta")
        MSAAligner.get_updated_alignment(locus_name=self.prg_builder.locus_name,
                                         previous_alignment=Path(previous_msa_filename),
                                         new_sequences=Path(new_sequences_filename),
                                         prefix=self.prg_builder.prefix)

        # update the alignment
        self.alignment = load_alignment_file(new_sequences_filename, "fasta")

        # update the derived attributes
        self._set_derived_helper_attributes()

        # update the recursion tree
        self._children = self._get_children()

        # reset the new sequences
        self.new_sequences = None
    ##################################################################################

    #####################################################################################################
    #  Clustering methods
    def _get_subalignments_by_clustering(self):
        id_lists = kmeans_cluster_seqs_in_interval(
            [self.interval.start, self.interval.stop],
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
        if self.interval is not None:
            sub_alignment = sub_alignment[:, self.interval.start : self.interval.stop + 1]
        return sub_alignment
    #####################################################################################################




class PrgBuilderSingleClusterNode(PrgBuilderRecursiveTreeNode):
    def __init__(self,
                 nesting_level,
                 alignment,
                 interval,
                 parent,
                 prg_builder):
        super().__init__(
            nesting_level,
            alignment,
            interval,
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
        logging.info(
            "match intervals: %s; non_match intervals: %s",
            self.match_intervals,
            self.non_match_intervals,
        )

    def _get_children(self):
        # stop conditions
        single_match_interval = (len(self.all_intervals) == 1) and (self.all_intervals[0] in self.match_intervals)
        max_nesting_level_reached = self.nesting_level == self.prg_builder.max_nesting
        small_variant_site = self.interval.stop - self.interval.start <= self.prg_builder.min_match_length
        if single_match_interval or max_nesting_level_reached or small_variant_site:
            return list()

        children = []
        for interval in self.all_intervals:
            if interval in self.match_intervals:
                # all seqs are not necessarily exactly the same: some can have 'N'
                # thus still process all of them, to get the one with no 'N'.
                sub_alignment = self.alignment[:, interval.start : interval.stop + 1]
                seqs = get_interval_seqs(sub_alignment)
                assert len(seqs) == 1, "Got >1 filtered sequences in match interval"
                child = PrgBuilderSingleClusterNode(
                    nesting_level=self.nesting_level + 1,
                    alignment=sub_alignment,
                    interval=interval,
                    parent=self,
                    prg_builder=self.prg_builder
                )
                children.append(child)
            else:
                child = PrgBuilderMultiClusterNode(
                    nesting_level=self.nesting_level + 1,
                    alignment=self.alignment,
                    interval=interval,
                    parent=self,
                    prg_builder=self.prg_builder
                )
                children.append(child)

        return children

    ##################################################################################
    # traversal methods
    # TODO memoize this?
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
                or column == {"-"}
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
                self.prg_builder.update_leaves_index(interval=(start_index, end_index), node=self)
            else:
                # TODO: disallow this for now?
                # TODO: can't it works for small variant site and nested variants
                # Add the variant seqs to the prg.
                site_num = self.prg_builder.get_next_site_num()
                prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")
                for seq_index, seq in enumerate(seqs):
                    site_num_for_this_seq = (site_num + 1) if (seq_index < len(seqs) - 1) else site_num
                    prg_as_list.extend(f"{seq}{delim_char}{site_num_for_this_seq}{delim_char}")

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
    # update methods
    # TODO: ignore for now
    def add_seqs_to_batch_update(self, new_sequences: List[Seq]):
        first_multi_cluster_parent = self._get_first_multi_cluster_parent()
        no_multi_cluster_parent_found = first_multi_cluster_parent is None
        if no_multi_cluster_parent_found:
            # update is either on leaf or root, let's avoid update on root, and update this leaf
            if self.new_sequences is None:
                self.new_sequences = []
            self.new_sequences.extend(new_sequences)
            self.new_sequences = list(set(self.new_sequences))  # dedup
        else:
            first_multi_cluster_parent.add_seqs_to_batch_update(new_sequences)

    def batch_update(self):
        no_update_to_be_done = self.new_sequences is None or len(self.new_sequences) == 0
        if no_update_to_be_done:
            return

        first_multi_cluster_parent = self._get_first_multi_cluster_parent()
        no_multi_cluster_parent_found = first_multi_cluster_parent is None
        if no_multi_cluster_parent_found:
            # update is either on leaf or root, let's avoid update on root, and update this leaf
            self._update_leaf()
        else:
            first_multi_cluster_parent.batch_update()


    def _update_leaf(self):
        # update the MSA
        previous_msa_filename = f"{self.prg_builder.prefix}.previous_msa.fa"
        with open(previous_msa_filename, "w") as previous_msa_handler:
            SeqIO.write(self.alignment, previous_msa_handler, "fasta")
        new_sequences_filename = f"{self.prg_builder.prefix}.new_sequences.fa"
        with open(new_sequences_filename, "w") as new_sequences_handler:
            SeqIO.write(self.new_sequences, new_sequences_handler, "fasta")
        MSAAligner.get_updated_alignment(locus_name=self.prg_builder.locus_name,
                                         previous_alignment=Path(previous_msa_filename),
                                         new_sequences=Path(new_sequences_filename),
                                         prefix=self.prg_builder.prefix)

        # update the alignment
        self.alignment = load_alignment_file(new_sequences_filename, "fasta")

        # update the other fields
        self._set_derived_helper_attributes()

        # reset the new sequences
        self.new_sequences = None

    def _get_first_multi_cluster_parent(self):
        node = self
        while True:
            no_multi_cluster_parent_found = node is None
            if no_multi_cluster_parent_found:
                return None
            is_a_multi_cluster_parent = isinstance(node, PrgBuilderMultiClusterNode)
            if is_a_multi_cluster_parent:
                return node
            node = node.parent
    ##################################################################################

class PrgBuilder(object):
    """
    Prg builder based from a multiple sequence alignment.
    Note min_match_length must be strictly greater than max_nesting + 1.
    """

    def __init__(
        self,
        locus_name,
        prefix,
        msa_file,
        alignment_format,
        max_nesting,
        min_match_length
    ):
        self.locus_name = locus_name
        self.prefix = prefix
        self.msa_file = msa_file
        self.alignment_format = alignment_format
        self.max_nesting = max_nesting
        self.min_match_length = min_match_length
        self.leaves_index = {}

        alignment = load_alignment_file(msa_file, alignment_format)
        root_interval = Interval(IntervalType.Root, 0, alignment.get_alignment_length() - 1)
        self._root = PrgBuilderSingleClusterNode(nesting_level=1,
                                                 alignment=alignment,
                                                 interval=root_interval,
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

    def update_leaves_index(self, interval, node):
        self.leaves_index[interval] = node

    def serialize(self, filename):
        with open(filename, "wb") as filehandler:
            pickle.dump(self, filehandler)

    @classmethod
    def deserialize(self, filename) -> "PrgBuilder":
        with open(filename, "rb") as filehandler:
            return pickle.load(filehandler)
