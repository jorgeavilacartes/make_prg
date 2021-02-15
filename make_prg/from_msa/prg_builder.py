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



class PrgBuilderMultiClusterNode(object):
    def __init__(self,
                 nesting_level,
                 alignments,
                 interval,
                 prg_builder):
        self.nesting_level = nesting_level
        self.alignments = alignments
        self.interval = interval
        self.prg_builder = prg_builder
        self._children = self._get_children()

    def _get_children(self):
        children = []
        for alignment in self.alignments:
            child = PrgBuilderSingleClusterNode(
                nesting_level=self.nesting_level,
                alignment=alignment,
                interval=self.interval,
                prg_builder=self.prg_builder
            )
            children.append(child)
        return children

    def _preorder_traversal_to_build_prg(self, delim_char) -> str:
        site_num = self.prg_builder.get_next_site_num()
        prg = f"{delim_char}{site_num}{delim_char}"
        
        for child_index, child in enumerate(self._children):
            site_num_to_separate_alleles = (site_num + 1) if (child_index < len(self._children) - 1) else site_num
            prg += child._preorder_traversal_to_build_prg(delim_char)
            prg += f"{delim_char}{site_num_to_separate_alleles}{delim_char}"
        return prg


class PrgBuilderSingleClusterNode(object):
    def __init__(self,
                 nesting_level,
                 alignment,
                 interval,
                 prg_builder):
        self.nesting_level = nesting_level
        self.alignment = alignment
        self.interval = interval
        self.prg_builder = prg_builder

        self.consensus = self.get_consensus(self.alignment)
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

        # generate recursion tree
        self._children = self._get_children()

    def _preorder_traversal_to_build_prg(self, delim_char) -> str:
        is_leaf_node = len(self._children) == 0
        if is_leaf_node:
            return self._get_prg(delim_char)
        else:
            prg = ""
            for child in self._children:
                child_prg = child._preorder_traversal_to_build_prg(delim_char)
                prg += child_prg
            return prg

    @classmethod
    def get_consensus(cls, alignment: MSA):
        """ Produces a 'consensus string' from an MSA: at each position of the
        MSA, the string has a base if all aligned sequences agree, and a "*" if not.
        IUPAC ambiguous bases result in non-consensus and are later expanded in the prg.
        N results in consensus at that position unless they are all N."""
        consensus_string = ""
        for i in range(alignment.get_alignment_length()):
            column = set([record.seq[i] for record in alignment])
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

    @classmethod
    def get_sub_alignment_by_list_id(
        self, id_list: List[str], alignment: MSA, interval=None
    ):
        list_records = [record for record in alignment if record.id in id_list]
        sub_alignment = MSA(list_records)
        if interval is not None:
            sub_alignment = sub_alignment[:, interval[0] : interval[1] + 1]
        return sub_alignment

    def _get_children(self):
        """
        Get and add the children of this node in the recursion tree
        """

        # stop conditions
        single_match_interval = (len(self.all_intervals)==1) and (self.all_intervals[0] in self.match_intervals)
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
                    prg_builder=self.prg_builder
                )
                children.append(child)

            else:
                # cluster
                id_lists = kmeans_cluster_seqs_in_interval(
                    [interval.start, interval.stop],
                    self.alignment,
                    self.prg_builder.min_match_length,
                )
                list_sub_alignments = [
                    self.get_sub_alignment_by_list_id(
                        id_list, self.alignment, [interval.start, interval.stop]
                    )
                    for id_list in id_lists
                ]

                child = PrgBuilderMultiClusterNode(
                    nesting_level=self.nesting_level + 1,
                    alignments=list_sub_alignments,
                    interval=interval,
                    prg_builder=self.prg_builder
                )
                children.append(child)

        return children

    def _get_prg(self, delim_char):
        prg = ""
        for interval in self.all_intervals:
            sub_alignment = self.alignment[:, interval.start : interval.stop + 1]
            seqs = get_interval_seqs(sub_alignment)

            single_seq = len(seqs) == 1
            if single_seq:
                prg = seqs[0]
            else:
                # Add the variant seqs to the prg.
                site_num = self.prg_builder.get_next_site_num()
                prg += f"{delim_char}{site_num}{delim_char}"
                for seq_index, seq in enumerate(seqs):
                    site_num_for_this_seq = (site_num + 1) if (seq_index < len(seqs) - 1) else site_num
                    prg += seq
                    prg += f"{delim_char}{site_num_for_this_seq}{delim_char}"
        return prg

    @property
    def max_nesting_level_reached(self):
        max_nesting = []
        if self.subAlignedSeqs == {}:
            logging.debug(
                "self.subAlignedSeqs == {} at nesting level %d for interval %s",
                self.nesting_level,
                self.interval,
            )
            max_nesting.append(self.nesting_level)
        else:
            logging.debug(
                "self.subAlignedSeqs.keys(): %s", list(self.subAlignedSeqs.keys())
            )
            logging.debug(
                "self.subAlignedSeqs[self.subAlignedSeqs.keys()[0]]: %s",
                self.subAlignedSeqs[list(self.subAlignedSeqs.keys())[0]],
            )
            for interval_start in list(self.subAlignedSeqs.keys()):
                logging.debug("interval start: %d", interval_start)
                for subaseq in self.subAlignedSeqs[interval_start]:
                    logging.debug(
                        "type of subAlignedSeqs object in list: %s", type(subaseq)
                    )
                    recur = subaseq.max_nesting_level_reached
                    logging.debug(
                        "recur max level nesting returned: %d, which has type %s",
                        recur,
                        type(recur),
                    )
                    max_nesting.append(recur)
        m = max(max_nesting)
        logging.debug("found the max of %s is %d", max_nesting, m)
        return m

    @property
    def prop_in_match_intervals(self):
        length_match_intervals = 0
        for interval in self.match_intervals:
            length_match_intervals += interval.stop - interval.start + 1
        return length_match_intervals / float(self.length)

    @property
    def num_seqs(self):
        return len(self.alignment)






class PrgBuilder(object):
    """
    Prg builder based from a multiple sequence alignment.
    Note min_match_length must be strictly greater than max_nesting + 1.
    """

    def __init__(
        self,
        msa_file,
        alignment_format="fasta",
        max_nesting=2,
        min_match_length=3
    ):
        self.msa_file = msa_file
        self.alignment_format = alignment_format
        self.max_nesting = max_nesting
        self.min_match_length = min_match_length

        alignment = load_alignment_file(msa_file, alignment_format)
        root_interval = Interval(IntervalType.Root, 0, alignment.get_alignment_length() - 1)
        self._root = PrgBuilderSingleClusterNode(nesting_level=1, alignment=alignment, interval=root_interval,
                                                 prg_builder=self)


    def build_prg(self):
        self._site_num = 5
        self.prg = self._root._preorder_traversal_to_build_prg(delim_char=" ")

    def get_next_site_num(self):
        previous_site_num = self._site_num
        self._site_num+=2
        return previous_site_num