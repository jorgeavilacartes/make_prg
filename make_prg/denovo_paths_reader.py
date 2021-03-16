from typing import List, Tuple, Optional
from collections import defaultdict, deque
import re
from intervaltree.intervaltree import IntervalTree
import logging
import sys
from Bio import pairwise2
from prg_builder import *


class MLPathNode:
    def __init__(self, key, sequence):
        assert len(sequence) > 0, "Error: pandora produced a Local Node with empty sequence"
        self.key = key
        self.sequence = sequence
        self.start_index_in_linear_path = None  # this is set by class MLPath if needed
        self.end_index_in_linear_path = None    # this is set by class MLPath if needed

    def __str__(self):
        return f"{self.key} {self.sequence} {self.start_index_in_linear_path} {self.end_index_in_linear_path}"

    def __repr__(self):
        return str(self)


class DenovoVariant:
    def __init__(self, start_index_in_linear_path, ref, alt):
        assert ref != alt, "Error: pandora produced a ref == alt"
        assert start_index_in_linear_path >= 0, "Error: Pandora produced a negative index for variant pos"
        self.start_index_in_linear_path = start_index_in_linear_path
        self.end_index_in_linear_path = start_index_in_linear_path + len(ref)
        self.ref = ref
        self.alt = alt

    def get_mutated_sequence(self, node: MLPathNode):
        start_index_inside_node_sequence = self.start_index_in_linear_path-node.start_index_in_linear_path
        end_index_inside_node_sequence = start_index_inside_node_sequence + len(self.ref)

        ref_wrt_indexes = node.sequence[start_index_inside_node_sequence:end_index_inside_node_sequence]
        ref_is_consistent = self.ref == ref_wrt_indexes
        assert ref_is_consistent, f"Ref is not consistent for {self}. Node = {node}. ref_wrt_indexes = {ref_wrt_indexes}"

        mutated_sequence = node.sequence[:start_index_inside_node_sequence] + \
                           self.alt + node.sequence[end_index_inside_node_sequence:]
        return mutated_sequence


    @staticmethod
    def align(ref, alt):
        alignment = pairwise2.align.globalms(ref, alt, MATCH_SCORE, MISMATCH_SCORE, GAP_OPEN_SCORE, GAP_EXTEND_SCORE)
        empty_alignment = len(alignment) == 0
        if empty_alignment:
            #  this usually happens if ref or alt are empty, let's check
            ref_is_empty = len(ref) == 0
            alt_is_empty = len(alt) == 0

            # only one should be empty
            both_are_empty = ref_is_empty and alt_is_empty
            assert not both_are_empty
            both_are_not_empty = (not ref_is_empty) and (not alt_is_empty)
            assert not both_are_not_empty

            if ref_is_empty:
                return "-" * len(alt), alt
            elif alt_is_empty:
                return ref, "-" * len(ref)
            else:
                assert True, "Unreachable code"  # just to be sure
        else:
            alignment_is_unique = len(alignment) == 1
            assert alignment_is_unique
            alignment = alignment[0]
            return alignment.seqA, alignment.seqB

    def split_variant(self, ml_path_nodes: List[MLPathNode], ml_path_nodes_to_nb_of_bases) -> Optional[List["DenovoVariant"]]:
        """
        Split this variant into a list of variants WRT to how it is distributed along the ML path
        """
        # 1. align ref to alt
        alignment = self.align(self.ref, self.alt)
        ref_alignment = deque(alignment[0])
        alt_alignment = deque(alignment[1])

        # 2. split variants at the boundaries of the alignments
        split_variants = []
        current_index_in_linear_path = self.start_index_in_linear_path
        for ml_path_node in ml_path_nodes:
            sub_ref = []
            sub_alt = []
            assert ml_path_node in ml_path_nodes_to_nb_of_bases
            nb_of_bases_to_consume = ml_path_nodes_to_nb_of_bases[ml_path_node]
            current_start_in_linear_path = current_index_in_linear_path

            while nb_of_bases_to_consume > 0:
                ref_base = ref_alignment.popleft()
                if ref_base != "-":
                    sub_ref.append(ref_base)
                    current_index_in_linear_path += 1
                    nb_of_bases_to_consume -= 1

                alt_base = alt_alignment.popleft()
                if alt_base != "-":
                    sub_alt.append(alt_base)

            split_variant = DenovoVariant(current_start_in_linear_path,
                                          "".join(sub_ref),
                                          "".join(sub_alt))
            split_variants.append(split_variant)

        return split_variants

    def is_insertion_event(self):
        return len(self.ref) == 0

    def __str__(self):
        return f"{self.start_index_in_linear_path} {self.ref} {self.alt}"

    def __repr__(self):
        return str(self)


class MLPath:
    def __init__(self, ml_path_nodes: List[MLPathNode]):
        self.__ml_path_nodes = ml_path_nodes
        self.__ml_path_index = IntervalTree()
        start_index_in_linear_path = 0
        for ml_path_node_index, ml_path_node in enumerate(ml_path_nodes):
            end_index_in_linear_path_for_indexing = start_index_in_linear_path + len(ml_path_node.sequence)
            is_last_node = ml_path_node_index == len(ml_path_nodes)-1
            if is_last_node:
                end_index_in_linear_path_for_indexing += 1  # allows for insertions after the ML sequence

            self.__ml_path_index.addi(start_index_in_linear_path,
                                      end_index_in_linear_path_for_indexing,
                                      data=ml_path_node)
            ml_path_node.start_index_in_linear_path = start_index_in_linear_path
            ml_path_node.end_index_in_linear_path = start_index_in_linear_path + len(ml_path_node.sequence)
            start_index_in_linear_path += len(ml_path_node.sequence)

    def get_node_at_index(self, index):
        nodes = self.__ml_path_index[index]
        only_one_node_is_overlapped = len(nodes) == 1
        assert only_one_node_is_overlapped
        node = list(nodes)[0].data
        return node


class MLPathNodeWithVariantApplied:
    def __init__(self, ml_path_node, variant, mutated_node_sequence):
        assert ml_path_node.sequence != mutated_node_sequence
        self.key = ml_path_node.key
        self.ml_path_node = ml_path_node
        self.variant = variant
        self.mutated_node_sequence = mutated_node_sequence


class DenovoLocusInfo:
    def __init__(self, sample: str, locus: str, ml_path: MLPath, variants: List[DenovoVariant]):
        self.sample = sample
        self.locus = locus
        self.ml_path = ml_path
        self.variants = variants

    def get_nodes_with_variant_applied(self) -> List[MLPathNodeWithVariantApplied]:
        """
        Get ML path nodes with the variants applied.
        @return: all nodes with variant applied
        """
        nodes_with_variant_applied = []
        for variant in self.variants:
            if variant.is_insertion_event():
                # interval is empty
                ml_path_node = self.ml_path.get_node_at_index(variant.start_index_in_linear_path)
                ml_path_nodes = [ml_path_node]
                ml_path_nodes_to_nb_of_bases = {ml_path_node: 1}
            else:
                # we have a real interval
                ml_path_nodes = []
                ml_path_nodes_to_nb_of_bases = defaultdict(int)
                for index_in_linear_path in range(variant.start_index_in_linear_path, variant.end_index_in_linear_path):
                    ml_path_node = self.ml_path.get_node_at_index(index_in_linear_path)
                    if ml_path_node not in ml_path_nodes:
                        ml_path_nodes.append(ml_path_node)
                    ml_path_nodes_to_nb_of_bases[ml_path_node] += 1

            assert len(ml_path_nodes) > 0
            variant_goes_through_several_leaves = len(ml_path_nodes) > 1
            if variant_goes_through_several_leaves:
                print(f"Variant goes through several leaves:", file=sys.stderr)
                print(f"self.sample: {self.sample}", file=sys.stderr)
                print(f"self.locus: {self.locus}", file=sys.stderr)
                print(f"variant: {variant}", file=sys.stderr)
                print(f"ml_path_nodes: {ml_path_nodes}", file=sys.stderr)
                print(f"ml_path_nodes_to_nb_of_bases: {ml_path_nodes_to_nb_of_bases}", file=sys.stderr)
                split_variants = variant.split_variant(ml_path_nodes, ml_path_nodes_to_nb_of_bases)

                split_was_unsuccessful = split_variants is None
                if split_was_unsuccessful:
                    continue

                print(f"split_variants: {split_variants}", file=sys.stderr)
            else:
                split_variants = [variant]

            assert len(split_variants) == len(ml_path_nodes)

            for split_variant, ml_path_node in zip(split_variants, ml_path_nodes):
                node_with_mutated_variant = MLPathNodeWithVariantApplied(
                    ml_path_node=ml_path_node,
                    variant=split_variant,
                    mutated_node_sequence=split_variant.get_mutated_sequence(ml_path_node)
                )
                nodes_with_variant_applied.append(node_with_mutated_variant)

        return nodes_with_variant_applied


class DenovoPathsDB:
    def __read_ml_path(self, filehandler):
        line = filehandler.readline().strip()
        nb_of_nodes_in_ml_path = int(line.split()[0])
        ml_path = []
        for _ in range(nb_of_nodes_in_ml_path):
            line = filehandler.readline().strip()

            try:
                matches = self.__ml_path_regex.search(line)
                start_index = int(matches.group(1))
                end_index = int(matches.group(2)) + 1  # TODO: fix pandora to give us non-inclusive end intervals instead
                sequence = matches.group(3)
            except Exception as exc:
                print(f"Failed matching ML path regex to line: {line}")
                print(f"Exception: {str(exc)}")
                sys.exit(1)

            assert start_index < end_index
            no_sequence_node = len(sequence) == 0
            if no_sequence_node:
                continue

            ml_path.append(MLPathNode(
                key=(start_index, end_index),
                sequence=sequence))
        return MLPath(ml_path)

    def __read_variants(self, filehandler):
        line = filehandler.readline().strip()
        nb_of_variants = int(line.split()[0])
        variants = []
        for _ in range(nb_of_variants):
            line = filehandler.readline()[:-1]
            line_split = line.split('\t')

            start_index_in_linear_path = int(line_split[0]) - 1
            ref = line_split[1]
            alt = line_split[2]
            variants.append(DenovoVariant(start_index_in_linear_path=start_index_in_linear_path, ref=ref, alt=alt))
        return variants

    def __populate__locus_name_to_denovo_loci(self):
        # Example:
        # (0 [0, 110) ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCACCGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC)
        self.__ml_path_regex = re.compile("\(\d+ \[(\d+), (\d+)\) ([ACGT]*)\)")

        self.__locus_name_to_denovo_loci = defaultdict(list)
        with open(self.filename) as filehandler:
            # read nb_of_samples
            line = filehandler.readline().strip()
            nb_of_samples = int(line.split()[0])

            for sample_index in range(nb_of_samples):
                # read each sample
                line = filehandler.readline().strip()
                sample = line.split()[1]
                line = filehandler.readline().strip()
                nb_of_loci_in_sample = int(line.split()[0])
                for locus_index in range(nb_of_loci_in_sample):
                    # read each locus
                    locus = filehandler.readline().strip()
                    ml_path = self.__read_ml_path(filehandler)
                    variants = self.__read_variants(filehandler)
                    denovo_locus = DenovoLocusInfo(sample, locus, ml_path, variants)
                    self.__locus_name_to_denovo_loci[locus].append(denovo_locus)

    def __populate__locus_name_to_variant_nodes_with_mutation(self):
        self.__locus_name_to_variant_nodes_with_mutation = defaultdict(list)
        for locus_name, denovo_loci in self.__locus_name_to_denovo_loci.items():
            for denovo_locus in denovo_loci:
                nodes_with_variant_applied = denovo_locus.get_nodes_with_variant_applied()
                self.__locus_name_to_variant_nodes_with_mutation[locus_name].extend(nodes_with_variant_applied)

    @property
    def locus_name_to_variant_nodes_with_mutation(self):
        return self.__locus_name_to_variant_nodes_with_mutation

    def __init__(self, filename):
        self.filename = filename
        self.__populate__locus_name_to_denovo_loci()
        self.__populate__locus_name_to_variant_nodes_with_mutation()
