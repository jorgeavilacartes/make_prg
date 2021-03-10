from typing import List, Tuple
from collections import defaultdict
import re
from intervaltree.intervaltree import IntervalTree
import logging
import sys
import types


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

        # TODO: move this back to assert
        if not ref_is_consistent:
            raise RuntimeError("Ref is not consistent.\n"
                               f"self = {self}\n"
                               f"node = {node}\n"
                               f"ref_wrt_indexes = {ref_wrt_indexes}\n"
                               f"start_index_inside_node_sequence = {start_index_inside_node_sequence}\n"
                               f"end_index_inside_node_sequence = {end_index_inside_node_sequence}")
        # assert ref_is_consistent, f"Ref is not consistent for {self}. Node = {node}. ref_wrt_indexes = {ref_wrt_indexes}"

        mutated_sequence = node.sequence[:start_index_inside_node_sequence] + \
                           self.alt + node.sequence[end_index_inside_node_sequence:]
        return mutated_sequence

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

    def get_nodes_with_variant_applied(self) -> Tuple[List[MLPathNodeWithVariantApplied], List[DenovoVariant]]:
        """
        Get ML path nodes with the variants applied.
        @return: a pair. First element is all nodes with variant applied when ref goes through a single node.
        If ref goes through two or more nodes, we add the variant to the list returned in the last element,
        and ignore it for now
        """
        nodes_with_variant_applied_in_a_single_node = []
        variants_in_two_or_more_nodes = []
        for variant in self.variants:
            node_of_first_base = self.ml_path.get_node_at_index(variant.start_index_in_linear_path)
            if not variant.is_insertion_event():
                node_of_last_base = self.ml_path.get_node_at_index(variant.end_index_in_linear_path-1)
            else:
                node_of_last_base = node_of_first_base
            variant_in_a_single_node = node_of_first_base == node_of_last_base
            if variant_in_a_single_node:
                try:
                    node_with_mutated_variant = MLPathNodeWithVariantApplied(
                        ml_path_node=node_of_first_base,
                        variant=variant,
                        mutated_node_sequence=variant.get_mutated_sequence(node_of_first_base)
                    )
                    nodes_with_variant_applied_in_a_single_node.append(node_with_mutated_variant)
                except RuntimeError as exc:
                    print("Failed applying variants to node sequence", file=sys.stderr)
                    print(exc, file=sys.stderr)
            else:
                variants_in_two_or_more_nodes.append(variant)

        return nodes_with_variant_applied_in_a_single_node, variants_in_two_or_more_nodes


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
        number_of_valid_variants = 0
        self.__ignored_variants_due_to_spanning_multiple_nodes = []
        for locus_name, denovo_loci in self.__locus_name_to_denovo_loci.items():
            for denovo_locus in denovo_loci:
                nodes_with_variant_applied_in_a_single_node, variants_in_two_or_more_nodes = \
                        denovo_locus.get_nodes_with_variant_applied()
                self.__locus_name_to_variant_nodes_with_mutation[locus_name].extend(
                    nodes_with_variant_applied_in_a_single_node)
                number_of_valid_variants += len(nodes_with_variant_applied_in_a_single_node)
                self.__ignored_variants_due_to_spanning_multiple_nodes.extend(variants_in_two_or_more_nodes)
        logging.warning(f"Number of valid variants (spanning a single node): {number_of_valid_variants}")
        logging.warning(f"Number of variants ignored due to spanning multiple nodes: "
                     f"{len(self.__ignored_variants_due_to_spanning_multiple_nodes)}")
        ignored_variants_as_str = "\n".join([str(x) for x in self.__ignored_variants_due_to_spanning_multiple_nodes])
        logging.warning(f"Variants spanning multiple nodes: {ignored_variants_as_str}")

    @property
    def locus_name_to_variant_nodes_with_mutation(self):
        return self.__locus_name_to_variant_nodes_with_mutation

    def __init__(self, filename):
        self.filename = filename
        self.__populate__locus_name_to_denovo_loci()
        self.__populate__locus_name_to_variant_nodes_with_mutation()
