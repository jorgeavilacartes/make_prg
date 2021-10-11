from typing import List, Set, Optional
from loguru import logger

from make_prg.from_msa import MSA
from make_prg.from_msa.cluster_sequences import kmeans_cluster_seqs_in_interval
from make_prg.seq_utils import (
    ambiguous_bases,
    remove_duplicates,  # TODO: quantify and remove dups
    get_interval_seqs,
    NONMATCH,
)
from make_prg.from_msa.interval_partition import IntervalPartitioner
from make_prg.prg_builder import PrgBuilder
from abc import ABC, abstractmethod
import copy
import numpy as np
from Bio.Seq import Seq
from Bio.AlignIO import MultipleSeqAlignment


class RecursiveTreeNode(ABC):
    def __init__(self, nesting_level: int, alignment: MultipleSeqAlignment, parent: Optional["RecursiveTreeNode"],
                 prg_builder: PrgBuilder):
        # set the basic attributes
        self.nesting_level: int = nesting_level
        self.parent: "RecursiveTreeNode" = parent
        self.prg_builder: PrgBuilder = prg_builder

        self.new_sequences: Set[str] = set()
        self.id: int = self.prg_builder.get_next_node_id()
        self.alignment: MultipleSeqAlignment = self._remove_gaps(alignment)
        self._set_derived_helper_attributes()

        # generate recursion tree
        self._children: List["RecursiveTreeNode"] = self._get_children()

    @staticmethod
    def _remove_gaps(alignment: MultipleSeqAlignment) -> MultipleSeqAlignment:
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

        gapless_alignment = MultipleSeqAlignment(gapless_records)
        return gapless_alignment

    @abstractmethod
    def _set_derived_helper_attributes(self):
        pass

    @abstractmethod
    def _get_children(self) -> List["RecursiveTreeNode"]:
        pass

    @abstractmethod
    def preorder_traversal_to_build_prg(self, prg_as_list: List[str], delim_char=" "):
        pass

    @abstractmethod
    def batch_update(self):
        pass


class MultiClusterNode(RecursiveTreeNode):
    def _set_derived_helper_attributes(self):
        pass  # nothing to set here

    def _get_children(self) -> List["RecursiveTreeNode"]:
        # each child is a PrgBuilderSingleClusterNode for each cluster subalignment
        cluster_subalignments = self._get_subalignments_by_clustering()
        children = []
        for alignment in cluster_subalignments:
            child = SingleClusterNode(
                nesting_level=self.nesting_level,
                alignment=alignment,
                parent=self,
                prg_builder=self.prg_builder
            )
            children.append(child)
        return children

    ##################################################################################
    # traversal methods
    def preorder_traversal_to_build_prg(self, prg_as_list: List[str], delim_char=" "):
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

    def batch_update(self):
        assert False, "Update was called on a non-leaf node"  # is an assertion because this a dev error


class SingleClusterNode(RecursiveTreeNode):
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

    def _get_children(self) -> List["RecursiveTreeNode"]:
        # base cases / stop conditions
        single_match_interval = (len(self.all_intervals) == 1) and (
            self.all_intervals[0] in self.match_intervals
        )
        max_nesting_level_reached = self.nesting_level == self.prg_builder.max_nesting
        small_variant_site = self.alignment.get_alignment_length() < self.prg_builder.min_match_length

        if single_match_interval or max_nesting_level_reached or small_variant_site:
            return list()

        children = []
        for interval in self.all_intervals:
            sub_alignment = self.alignment[:, interval.start : interval.stop + 1]
            is_a_match_interval = interval in self.match_intervals
            if is_a_match_interval:
                # all seqs are not necessarily exactly the same: some can have 'N'
                # thus still process all of them, to get the one with no 'N'.
                seqs = get_interval_seqs(sub_alignment)
                assert len(seqs) == 1, "Got >1 filtered sequences in match interval"
                subclass = SingleClusterNode
            else:
                subclass = MultiClusterNode

            child = subclass(
                nesting_level=self.nesting_level + 1,
                alignment=sub_alignment,
                parent=self,
                prg_builder=self.prg_builder
            )
            children.append(child)

        return children

    ##################################################################################
    # traversal methods
    def preorder_traversal_to_build_prg(self, prg_as_list: List["str"], delim_char=" "):
        is_leaf_node = len(self._children) == 0
        if is_leaf_node:
            self._get_prg(prg_as_list, delim_char)
        else:
            for child in self._children:
                child.preorder_traversal_to_build_prg(prg_as_list, delim_char)
    ##################################################################################

    ##################################################################################
    # helpers
    def _get_consensus(self) -> str:
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

    def _get_prg(self, prg_as_list: List[str], delim_char=" "):
        for interval in self.all_intervals:
            sub_alignment = self.alignment[:, interval.start:interval.stop + 1]
            seqs = get_interval_seqs(sub_alignment)

            single_seq = len(seqs) == 1
            if single_seq:
                start_index = len(prg_as_list)
                prg_as_list.extend(seqs[0])
                end_index = len(prg_as_list) + 1
                self.prg_builder.update_leaves_index(start_index, end_index, node=self)
            else:
                # Add the variant seqs to the prg
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
                    prg_as_list.extend(
                        f"{delim_char}{site_num_for_this_seq}{delim_char}"
                    )
    ##################################################################################

    ##################################################################################
    # update methods
    def add_seq_to_batch_update(self, new_sequence: str):
        self.new_sequences.add(new_sequence)

    def batch_update(self):
        no_update_to_be_done = len(self.new_sequences) == 0
        if no_update_to_be_done:
            return
        self._update_leaf()

    def _update_leaf(self):
        logger.debug(
            f"Updating MSA for {self.prg_builder.locus_name}, node {self.id}..."
        )

        self.alignment = self.prg_builder.aligner.get_updated_alignment(
            current_alignment=self.alignment,
            new_sequences=self.new_sequences
        )

        # update the other fields
        self._set_derived_helper_attributes()

        # reset the new sequences
        self.new_sequences = set()

        # regenerate recursion tree
        self._children = self._get_children()
    ##################################################################################
