from typing import List, Tuple

from loguru import logger

from make_prg.from_msa import MSA
from make_prg.io_utils import load_alignment_file
from make_prg.from_msa.cluster_sequences import kmeans_cluster_seqs_in_interval
from make_prg.seq_utils import (
    ambiguous_bases,
    remove_duplicates,  # TODO: quantify and remove dups
    get_interval_seqs,
    NONMATCH,
)
from make_prg.from_msa.interval_partition import (
    IntervalPartitioner,
    Interval,
    IntervalType,
)
import pickle
from pathlib import Path
import os
import shlex
import time
from Bio import SeqIO
from abc import ABC, abstractmethod
import uuid
import shutil
import subprocess
import copy
import numpy as np
from Bio.Seq import Seq
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict

MATCH_SCORE = 2
MISMATCH_SCORE = -1
GAP_OPEN_SCORE = -4
GAP_EXTEND_SCORE = -2


class MSAAligner(ABC):
    def __init__(self, executable: str):
        self.executable = executable

    @abstractmethod
    def is_executable(self) -> bool:
        pass

    @abstractmethod
    def get_updated_alignment(
        self,
        leaf_name: str,
        previous_alignment: Path,
        new_sequences: List[str],
        temp_prefix: Path,
    ) -> Path:
        pass


class MSAAlignerMAFFT(MSAAligner):
    def is_executable(self) -> bool:
        return shutil.which(self.executable) is not None

    def get_updated_alignment(
        self,
        leaf_name: str,
        previous_alignment: Path,
        new_sequences: List[str],
        temp_prefix: Path,
    ) -> Path:
        new_sequences_filename = temp_prefix / "new_sequences.fa"
        with open(new_sequences_filename, "w") as new_sequences_handler:
            for index_new_seq, new_seq in enumerate(new_sequences):
                print(
                    f">Denovo_path_{index_new_seq}_random_id_{uuid.uuid4()}",
                    file=new_sequences_handler,
                )
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
                self.executable,
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
            args,
            stderr=subprocess.PIPE,
            encoding="utf-8",
            shell=True,
            env=env,
        )
        exit_code = process.wait()
        shutil.rmtree(mafft_tmpdir)
        if exit_code != 0:
            raise RuntimeError(
                f"Failed to execute MAFFT for {leaf_name} due to the following error:\n"
                f"{process.stderr.read()}"
            )
        stop = time.time()
        runtime = stop - start
        logger.debug(f"MAFFT update runtime for {leaf_name} in seconds: {runtime:.3f}")

        return new_msa


class PrgBuilderRecursiveTreeNode(ABC):
    def __init__(self, nesting_level, alignment, parent, prg_builder, mafft: str):
        # set the basic attributes
        self.id = prg_builder.get_next_node_id()
        self.nesting_level = nesting_level
        self.parent = parent
        self.prg_builder = prg_builder
        self.new_sequences = None

        self.alignment = self.remove_gaps(alignment)
        self.mafft = mafft

        self._set_derived_helper_attributes()

        self.prg_builder.all_nodes.append(self)

        is_root = self.id == 0
        if not is_root:
            self.prg_builder.graph.add_node(self.id)

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
    def batch_update(self, temp_prefix, mafft: str = "mafft"):
        pass


class PrgBuilderMultiClusterNode(PrgBuilderRecursiveTreeNode):
    def __init__(self, nesting_level, alignment, parent, prg_builder, mafft):
        super().__init__(nesting_level, alignment, parent, prg_builder, mafft=mafft)

    def _set_derived_helper_attributes(self):
        pass  # nothing to set here

    def _get_children(self):
        # each child is a PrgBuilderSingleClusterNode for each cluster subalignment
        cluster_subalignments = self._get_subalignments_by_clustering()
        children = []
        for alignment in cluster_subalignments:
            child = PrgBuilderSingleClusterNode(
                nesting_level=self.nesting_level + 1,  # TODO FIXME: this was added just for drawing purposes, roll back
                alignment=alignment,
                parent=self,
                prg_builder=self.prg_builder,
                mafft=self.mafft,
            )
            children.append(child)
            self.prg_builder.graph.add_edge(self.id, child.id)
        return children

    ##################################################################################
    # traversal methods
    def preorder_traversal_to_build_prg(self, prg_as_list, delim_char):
        site_num = self.prg_builder.get_next_site_num()
        prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")

        for child_index, child in enumerate(self._children):
            site_num_to_separate_alleles = (
                (site_num + 1) if (child_index < len(self._children) - 1) else site_num
            )
            child.preorder_traversal_to_build_prg(prg_as_list, delim_char)
            prg_as_list.extend(
                f"{delim_char}{site_num_to_separate_alleles}{delim_char}"
            )

    ##################################################################################

    #####################################################################################################
    #  Clustering methods
    def _get_subalignments_by_clustering(self):
        id_lists = kmeans_cluster_seqs_in_interval(
            self.alignment,
            self.prg_builder.min_match_length,
        )
        list_sub_alignments = [
            self._get_sub_alignment_by_list_id(id_list) for id_list in id_lists
        ]
        return list_sub_alignments

    def _get_sub_alignment_by_list_id(self, id_list: List[str]):
        list_records = [record for record in self.alignment if record.id in id_list]
        sub_alignment = MSA(list_records)
        return sub_alignment

    #####################################################################################################

    def batch_update(self, temp_prefix, mafft: str = ""):
        assert False, "Update was called on a non-leaf node"


class PrgBuilderSingleClusterNode(PrgBuilderRecursiveTreeNode):
    def __init__(self, nesting_level, alignment, parent, prg_builder, mafft):
        super().__init__(nesting_level, alignment, parent, prg_builder, mafft=mafft)

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
        single_match_interval = (len(self.all_intervals) == 1) and (
            self.all_intervals[0] in self.match_intervals
        )
        max_nesting_level_reached = self.nesting_level == self.prg_builder.max_nesting
        small_variant_site = (
            self.alignment.get_alignment_length() < self.prg_builder.min_match_length
        )
        if single_match_interval or max_nesting_level_reached or small_variant_site:
            return list()

        children = []
        for interval in self.all_intervals:
            sub_alignment = self.alignment[:, interval.start : interval.stop + 1]

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
                prg_builder=self.prg_builder,
                mafft=self.mafft,
            )
            children.append(child)
            edge_from_root = self.id == 0
            if not edge_from_root:
                self.prg_builder.graph.add_edge(self.id, child.id)
            has_siblings = len(children) > 1
            if has_siblings:
                self.prg_builder.graph.add_edge(children[-2].id, child.id)

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
        """Produces a 'consensus string' from an MSA: at each position of the
        MSA, the string has a base if all aligned sequences agree, and a "*" if not.
        IUPAC ambiguous bases result in non-consensus and are later expanded in the prg.
        N results in consensus at that position unless they are all N."""
        consensus_string = ""
        for i in range(self.alignment.get_alignment_length()):
            column = set([record.seq[i] for record in self.alignment])
            column = column.difference({"N"})
            if len(ambiguous_bases.intersection(column)) > 0 or len(column) != 1:
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
                end_index = len(prg_as_list) + 1
                self.prg_builder.update_leaves_index(start_index, end_index, node=self)
                self.prg_builder.leaf_id_to_seq[self.id] = seqs[0]
            else:
                # Add the variant seqs to the prg.
                site_num = self.prg_builder.get_next_site_num()
                prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")
                for seq_index, seq in enumerate(seqs):
                    site_num_for_this_seq = (
                        (site_num + 1) if (seq_index < len(seqs) - 1) else site_num
                    )
                    start_index = len(prg_as_list)
                    prg_as_list.extend(seq)
                    end_index = len(prg_as_list) + 1
                    self.prg_builder.update_leaves_index(
                        start_index, end_index, node=self
                    )
                    self.prg_builder.leaf_id_to_seq[self.id] = seq
                    prg_as_list.extend(
                        f"{delim_char}{site_num_for_this_seq}{delim_char}"
                    )

    ##################################################################################

    ##################################################################################
    # update methods
    def add_seq_to_batch_update(self, new_sequence: str):
        if self.new_sequences is None:
            self.new_sequences = set()
        self.new_sequences.add(new_sequence)

    def batch_update(self, temp_prefix, mafft: str = "mafft"):
        self.mafft = mafft
        no_update_to_be_done = (
            self.new_sequences is None or len(self.new_sequences) == 0
        )
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

        logger.debug(
            f"Updating MSA for {self.prg_builder.locus_name}, node {self.id}..."
        )
        msa_aligner = MSAAlignerMAFFT(self.mafft)
        logger.debug("Detecting mafft...")
        mafft_is_runnable = msa_aligner.is_executable()

        if mafft_is_runnable:
            logger.debug("mafft detected!")
        else:
            raise RuntimeError(
                f"Could not execute mafft using `{msa_aligner.executable}`"
            )

        new_msa = msa_aligner.get_updated_alignment(
            leaf_name=f"{self.prg_builder.locus_name}, node {self.id}",
            previous_alignment=previous_msa_filename,
            new_sequences=list(self.new_sequences),
            temp_prefix=temp_prefix,
        )

        # update the alignment
        self.alignment = load_alignment_file(str(new_msa), "fasta")

        # update the other fields
        self._set_derived_helper_attributes()

        # reset the new sequences
        self.new_sequences = None

        # regenerate recursion tree
        self._children = self._get_children()

    ##################################################################################


class LeafNotFoundException(Exception):
    pass


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
        min_match_length,
        mafft: str = "",
    ):
        self.locus_name = locus_name
        self.msa_file = msa_file
        self.alignment_format = alignment_format
        self.max_nesting = max_nesting
        self.min_match_length = min_match_length
        self.leaves_index = {}
        self.node_id = 0
        self.mafft = mafft
        self.all_nodes = []
        self.graph = nx.DiGraph()
        self.leaf_id_to_seq = {}

        alignment = load_alignment_file(msa_file, alignment_format)
        self._root = PrgBuilderSingleClusterNode(
            nesting_level=1,
            alignment=alignment,
            parent=None,
            prg_builder=self,
            mafft=self.mafft,
        )

    def build_prg(self):
        self._site_num = 5
        prg_as_list = []
        self._root.preorder_traversal_to_build_prg(prg_as_list, delim_char=" ")
        prg = "".join(prg_as_list)
        return prg

    def get_next_site_num(self):
        previous_site_num = self._site_num
        self._site_num += 2
        return previous_site_num

    def get_next_node_id(self):
        self.node_id += 1
        return self.node_id - 1

    def update_leaves_index(self, start_index: int, end_index: int, node):
        interval = (start_index, end_index)
        self.leaves_index[interval] = node

    def get_node_given_interval(self, interval: Tuple[int, int]):
        # TODO: move this back to assert once is solved
        interval_is_indexed = interval in self.leaves_index
        if not interval_is_indexed:
            raise LeafNotFoundException(
                f"Queried interval {interval} does not exist in leaves index for locus {self.locus_name}"
            )

        # assert interval in self.leaves_index, \
        #     f"Fatal error: Queried interval {interval} does not exist in leaves index for locus {self.locus_name}"

        return self.leaves_index[interval]

    def serialize(self, filename):
        with open(filename, "wb") as filehandler:
            pickle.dump(self, filehandler)

    @staticmethod
    def deserialize(filename):
        with open(filename, "rb") as filehandler:
            return pickle.load(filehandler)

    def output_graph(self, filename):
        nesting_level_to_node_ids = defaultdict(list)
        for node in self.all_nodes:
            nesting_level_to_node_ids[node.nesting_level].append(node.id)


        plt.figure(figsize=(20, 10))
        a_graph = nx.drawing.nx_agraph.to_agraph(self.graph)
        for node in a_graph.nodes():
            node.attr["label"] = self.leaf_id_to_seq.get(int(node.name), node.name)

        for nesting_level, nodes in nesting_level_to_node_ids.items():
            a_graph.add_subgraph(nodes, rank="same")
        a_graph.layout(prog="dot")
        a_graph.draw(filename)

        a_graph.draw(f"{filename}.dot")

        ## See original saving of each node location to node_pos
        # return node_pos
        # pos = nx.drawing.nx_agraph.graphviz_layout(self.graph, prog='dot', args='-Gnodesep=2 -layout=layout_reingold_tilford')
        # nx.draw(self.graph, pos=pos, with_labels=True, font_weight='bold', arrows=True)
        # plt.savefig(filename)


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

    def to_absolute_paths(self, parent: Path):
        for locus_name, pickle_file in self.locus_name_to_pickle_files.items():
            absolute_filepath = parent / pickle_file
            self.locus_name_to_pickle_files[locus_name] = str(absolute_filepath)
