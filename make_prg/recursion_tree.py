from typing import List, Set, Optional
from loguru import logger
from make_prg.from_msa import MSA
from make_prg.from_msa.cluster_sequences import kmeans_cluster_seqs
from make_prg.utils.seq_utils import (
    SequenceExpander,
    remove_columns_full_of_gaps_from_MSA,
    get_consensus_from_MSA,
    get_number_of_unique_ungapped_sequences,
    get_number_of_unique_gapped_sequences
)
from make_prg.from_msa.interval_partition import IntervalPartitioner, Interval
from abc import ABC, abstractmethod


class RecursiveTreeNode(ABC):
    def __init__(self, nesting_level: int, alignment: MSA, parent: Optional["RecursiveTreeNode"],
                 prg_builder: "PrgBuilder", force_no_child: bool = False):
        self.nesting_level: int = nesting_level
        self.alignment: MSA = remove_columns_full_of_gaps_from_MSA(alignment)
        self.parent: "RecursiveTreeNode" = parent
        self.prg_builder: "PrgBuilder" = prg_builder
        self.force_no_child = force_no_child
        self.id: int = self.prg_builder.get_next_node_id()
        self.new_sequences: Set[str] = set()

        # generate recursion tree
        self._init_pre_recursion_attributes()
        self._children: List["RecursiveTreeNode"] = self._get_children()

    @property
    def children(self):
        return self._children

    @abstractmethod
    def _init_pre_recursion_attributes(self):
        pass

    @abstractmethod
    def _get_children(self) -> List["RecursiveTreeNode"]:
        pass

    @abstractmethod
    def preorder_traversal_to_build_prg(self, prg_as_list: List[str], delim_char: str = " "):
        pass

    def is_leaf(self) -> bool:
        return len(self.children) == 0


class MultiClusterNode(RecursiveTreeNode):
    def __init__(self, nesting_level: int, alignment: MSA, parent: Optional["RecursiveTreeNode"],
                 prg_builder: "PrgBuilder", force_no_child: bool = False):
        super().__init__(nesting_level, alignment, parent, prg_builder, force_no_child)
        assert not self.is_leaf(), "Multicluster nodes should never be leaves"

    def _init_pre_recursion_attributes(self):
        pass  # nothing to init here

    def _get_children(self) -> List["RecursiveTreeNode"]:
        # each child is a PrgBuilderSingleClusterNode for each cluster subalignment
        cluster_subalignments = self._get_subalignments_by_clustering()
        no_clustering_was_done = len(cluster_subalignments) == 1
        children = []
        for alignment in cluster_subalignments:
            child = SingleClusterNode(
                nesting_level=self.nesting_level,
                alignment=alignment,
                parent=self,
                prg_builder=self.prg_builder,
                force_no_child=no_clustering_was_done
            )
            children.append(child)
        return children

    ##################################################################################
    # traversal methods
    def preorder_traversal_to_build_prg(self, prg_as_list: List[str], delim_char: str = " "):
        site_num = self.prg_builder.get_next_site_num()
        prg_as_list.extend(f"{delim_char}{site_num}{delim_char}")

        for child_index, child in enumerate(self.children):
            site_num_to_separate_alleles = (
                (site_num + 1) if (child_index < len(self.children) - 1) else site_num
            )
            child.preorder_traversal_to_build_prg(prg_as_list, delim_char)
            prg_as_list.extend(
                f"{delim_char}{site_num_to_separate_alleles}{delim_char}"
            )
    ##################################################################################

    #####################################################################################################
    #  clustering methods
    def _get_subalignments_by_clustering(self) -> List[MSA]:
        clustered_ids = kmeans_cluster_seqs(
            self.alignment,
            self.prg_builder.min_match_length
        )
        list_sub_alignments = [
            self._get_sub_alignment_by_list_id(clustered_id) for clustered_id in clustered_ids
        ]
        return list_sub_alignments

    def _get_sub_alignment_by_list_id(self, id_list: List[str]) -> MSA:
        list_records = [record for record in self.alignment if record.id in id_list]
        sub_alignment = MSA(list_records)
        return sub_alignment
    #####################################################################################################


class SingleClusterNode(RecursiveTreeNode):
    def _init_pre_recursion_attributes(self):
        self.consensus: str = get_consensus_from_MSA(self.alignment)
        self.length: int = len(self.consensus)
        (
            self.match_intervals,
            self.non_match_intervals,
            self.all_intervals,
        ) = IntervalPartitioner(
            self.consensus, self.prg_builder.min_match_length, self.alignment
        ).get_intervals()

    def _infer_if_this_node_has_no_child(self, alignment):
        if self.force_no_child:
            return True

        single_match_interval = (len(self.all_intervals) == 1) and (
                self.all_intervals[0] in self.match_intervals
        )
        if single_match_interval:
            return True

        max_nesting_level_reached = self.nesting_level == self.prg_builder.max_nesting
        if max_nesting_level_reached:
            return True

        small_variant_site = self.alignment.get_alignment_length() < self.prg_builder.min_match_length
        if small_variant_site:
            return True

        num_unique_nongapped_seqs = get_number_of_unique_ungapped_sequences(alignment)
        too_few_unique_sequences = num_unique_nongapped_seqs <= 2
        if too_few_unique_sequences:
            return True

        num_unique_gapped_seqs = get_number_of_unique_gapped_sequences(alignment)
        assert num_unique_nongapped_seqs <= num_unique_gapped_seqs
        alignment_has_ambiguity = num_unique_nongapped_seqs < num_unique_gapped_seqs
        if alignment_has_ambiguity:
            return True

        return False

    def _infer_if_should_not_cluster(self, interval: Interval, alignment: MSA):
        is_a_match_interval = interval in self.match_intervals
        if is_a_match_interval:
            return True

        num_unique_nongapped_seqs = get_number_of_unique_ungapped_sequences(alignment)
        too_few_unique_sequences = num_unique_nongapped_seqs <= 2
        if too_few_unique_sequences:
            return True

        num_unique_gapped_seqs = get_number_of_unique_gapped_sequences(alignment)
        assert num_unique_nongapped_seqs <= num_unique_gapped_seqs
        alignment_has_ambiguity = num_unique_nongapped_seqs < num_unique_gapped_seqs
        if alignment_has_ambiguity:
            return True

        clustered_ids = kmeans_cluster_seqs(alignment, self.prg_builder.min_match_length)
        single_cluster_found = len(clustered_ids) == 1
        if single_cluster_found:
            return True

        return False

    def _get_children(self) -> List["RecursiveTreeNode"]:
        node_has_no_child = self._infer_if_this_node_has_no_child(self.alignment)
        if node_has_no_child:
            return list()

        children = []
        for interval in self.all_intervals:
            sub_alignment = self.alignment[:, interval.start : interval.stop + 1]
            sub_alignment_will_not_be_reclustered = self._infer_if_should_not_cluster(interval, sub_alignment)
            if sub_alignment_will_not_be_reclustered:
                subclass = SingleClusterNode
            else:
                subclass = MultiClusterNode
            child = subclass(
                nesting_level=self.nesting_level + 1,
                alignment=sub_alignment,
                parent=self,
                prg_builder=self.prg_builder,
                force_no_child=sub_alignment_will_not_be_reclustered
            )
            children.append(child)

        return children

    def _get_prg(self, prg_as_list: List[str], delim_char: str = " "):
        for interval in self.all_intervals:
            sub_alignment = self.alignment[:, interval.start:interval.stop + 1]
            seqs = SequenceExpander.get_expanded_sequences(sub_alignment)

            single_seq = len(seqs) == 1
            if single_seq:
                start_index = len(prg_as_list)
                prg_as_list.extend(seqs[0])
                end_index = len(prg_as_list)
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
                    end_index = len(prg_as_list)
                    self.prg_builder.update_leaves_index(
                        start_index, end_index, node=self
                    )
                    prg_as_list.extend(
                        f"{delim_char}{site_num_for_this_seq}{delim_char}"
                    )

    ##################################################################################
    # traversal methods
    def preorder_traversal_to_build_prg(self, prg_as_list: List["str"], delim_char: str = " "):
        if self.is_leaf():
            self._get_prg(prg_as_list, delim_char)
        else:
            for child in self.children:
                child.preorder_traversal_to_build_prg(prg_as_list, delim_char)
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

        an_aligner_was_given = self.prg_builder.aligner is not None
        # this is an assertion as it is the dev responsibility to ensure an aligner is given if updates are to be done
        assert an_aligner_was_given, "Cannot make updates without an Multiple Sequence Aligner."

        self.alignment = self.prg_builder.aligner.get_updated_alignment(
            current_alignment=self.alignment,
            new_sequences=self.new_sequences
        )

        # reset the new sequences
        self.new_sequences = set()

        # regenerate recursion tree
        self._init_pre_recursion_attributes()
        self._children = self._get_children()
    ##################################################################################
