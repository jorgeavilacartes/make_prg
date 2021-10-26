from typing import List, Deque, TextIO, Dict, Tuple
from collections import defaultdict, deque
import re
from loguru import logger
from collections import Counter
from make_prg.utils.seq_utils import align, GAP
from make_prg.update.MLPath import MLPathNode, MLPath, EmptyMLPathSequence
from pathlib import Path
from itertools import groupby


class DenovoError(Exception):
    pass


class DenovoVariant:
    def __init__(self, start_index_in_linear_path: int, ref: str, alt: str):
        DenovoVariant._param_checking(start_index_in_linear_path, ref, alt)
        self.start_index_in_linear_path: int = start_index_in_linear_path
        self.end_index_in_linear_path: int = start_index_in_linear_path + len(ref)
        self.ref: str = ref
        self.alt: str = alt

    @staticmethod
    def _param_checking(start_index_in_linear_path: int, ref: str, alt: str):
        DenovoVariant._check_sequence_is_composed_of_ACGT_only(ref)
        DenovoVariant._check_sequence_is_composed_of_ACGT_only(alt)
        not_a_variant = ref == alt
        if not_a_variant:
            raise DenovoError(f"Found a variant where ref ({ref}) equals alt ({alt}), this is not a variant")

        negative_index_for_variant_pos = start_index_in_linear_path < 0
        if negative_index_for_variant_pos:
            raise DenovoError(f"Found a negative index for variant pos ({start_index_in_linear_path})")

    @staticmethod
    def _check_sequence_is_composed_of_ACGT_only(seq: str):
        sequence_is_composed_of_ACGT_only = all([base in "ACGT" for base in seq])
        if not sequence_is_composed_of_ACGT_only:
            raise DenovoError(f"Found a non-ACGT seq ({seq}) in a denovo variant")

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.start_index_in_linear_path, self.end_index_in_linear_path, self.ref, self.alt) == \
                   (other.start_index_in_linear_path, other.end_index_in_linear_path, other.ref, other.alt)
        else:
            return False

    def get_mutated_sequence(self, node: MLPathNode) -> str:
        node_is_compatible_with_this_variant = \
            node.start_index_in_linear_path <= self.start_index_in_linear_path and \
            self.end_index_in_linear_path <= node.end_index_in_linear_path
        if not node_is_compatible_with_this_variant:
            raise DenovoError(f"Node {node} is not compatible with variant {self}")

        start_index_inside_node_sequence = (
            self.start_index_in_linear_path - node.start_index_in_linear_path
        )
        end_index_inside_node_sequence = start_index_inside_node_sequence + len(
            self.ref
        )
        ref_wrt_indexes = node.sequence[
            start_index_inside_node_sequence:end_index_inside_node_sequence
        ]
        ref_is_consistent = self.ref == ref_wrt_indexes
        if not ref_is_consistent:
            raise DenovoError(f"Ref is not consistent for {self}. Node = {node}. ref_wrt_indexes = {ref_wrt_indexes}")

        mutated_sequence = (
            node.sequence[:start_index_inside_node_sequence]
            + self.alt
            + node.sequence[end_index_inside_node_sequence:]
        )
        return mutated_sequence

    def _split_variant_at_boundary_alignment(self, ml_path_nodes_it_goes_through: List[MLPathNode],
                                             ref_alignment: Deque[str], alt_alignment: Deque[str]
                                             ) -> List["DenovoVariant"]:
        split_variants = []
        current_index_in_linear_path = self.start_index_in_linear_path
        ml_path_node_to_count = Counter(ml_path_nodes_it_goes_through)
        deduplicated_ml_path_nodes_it_goes_through = [x[0] for x in groupby(ml_path_nodes_it_goes_through)]
        for ml_path_node_index, ml_path_node in enumerate(deduplicated_ml_path_nodes_it_goes_through):
            sub_ref = []
            sub_alt = []
            nb_of_bases_to_consume = ml_path_node_to_count[ml_path_node]
            current_start_in_linear_path = current_index_in_linear_path

            while nb_of_bases_to_consume > 0:
                ref_base = ref_alignment.popleft()
                if ref_base != GAP:
                    sub_ref.append(ref_base)
                    current_index_in_linear_path += 1
                    nb_of_bases_to_consume -= 1

                alt_base = alt_alignment.popleft()
                if alt_base != GAP:
                    sub_alt.append(alt_base)

            is_last_node = ml_path_node_index == len(deduplicated_ml_path_nodes_it_goes_through) - 1
            there_are_remaining_alt_bases = is_last_node and nb_of_bases_to_consume==0 and len(alt_alignment)>0
            if there_are_remaining_alt_bases:
                sub_alt.extend([alt_base for alt_base in alt_alignment if alt_base != GAP])

            sub_ref_seq = "".join(sub_ref)
            sub_alt_seq = "".join(sub_alt)
            sub_ref_and_alt_are_different = sub_ref_seq != sub_alt_seq
            if sub_ref_and_alt_are_different:
                split_variant = DenovoVariant(
                    current_start_in_linear_path, sub_ref_seq, sub_alt_seq
                )
                split_variants.append(split_variant)

        return split_variants

    def split_variant(
        self, ml_path_nodes_it_goes_through: List[MLPathNode]
    ) -> List["DenovoVariant"]:
        """
        Split this variant into a list of sub-variants WRT how it is distributed along the ML path
        @param: ml_path_nodes_it_goes_through: a list of MLPathNode, where the i-th MLPathNode is the node the i-th base
        of the variant goes through
        """
        each_base_is_covered_by_one_node = len(self.ref) == len(ml_path_nodes_it_goes_through)
        assert each_base_is_covered_by_one_node, "We have bases uncovered by nodes in split_variant()"

        nb_of_distinct_ml_path_nodes = len(set(ml_path_nodes_it_goes_through))
        variant_goes_through_only_one_leaf = nb_of_distinct_ml_path_nodes == 1
        if variant_goes_through_only_one_leaf:
            split_variants = [self]
            return split_variants

        # here, variant goes through several leaves
        alignment = align(self.ref, self.alt)
        ref_alignment = deque(alignment[0])
        alt_alignment = deque(alignment[1])
        split_variants = self._split_variant_at_boundary_alignment(
            ml_path_nodes_it_goes_through, ref_alignment, alt_alignment
        )
        return split_variants

    def is_insertion_event(self) -> bool:
        return len(self.ref) == 0

    def __str__(self):
        return f"{self.start_index_in_linear_path} {self.ref} {self.alt}"

    def __repr__(self):
        return str(self)


class UpdateData:
    def __init__(self, ml_path_node_key: Tuple[int, int], new_node_sequence: str):
        self.ml_path_node_key: Tuple[int, int] = ml_path_node_key
        self.new_node_sequence: str = new_node_sequence


class DenovoLocusInfo:
    def __init__(
        self, sample: str, locus: str, ml_path: MLPath, variants: List[DenovoVariant]
    ):
        self.sample: str = sample
        self.locus: str = locus
        self.ml_path: MLPath = ml_path
        self.variants: List[DenovoVariant] = variants

    def _get_ml_path_nodes_spanning_variant(self, variant: DenovoVariant) -> List[MLPathNode]:
        if variant.is_insertion_event():
            # interval is empty
            ml_path_node = self.ml_path.get_node_at_position(variant.start_index_in_linear_path)
            ml_path_nodes = [ml_path_node]
        else:
            # we have a real interval
            ml_path_nodes = []
            for position_in_linear_path in range(
                    variant.start_index_in_linear_path, variant.end_index_in_linear_path
            ):
                ml_path_node = self.ml_path.get_node_at_position(position_in_linear_path)
                ml_path_nodes.append(ml_path_node)

        return ml_path_nodes

    def get_update_data(self) -> List[UpdateData]:
        """
        Get a list of updates to be done to ML path nodes
        """
        update_data_list = []
        for variant in self.variants:
            ml_path_nodes = self._get_ml_path_nodes_spanning_variant(variant)
            split_variants = variant.split_variant(ml_path_nodes)

            for split_variant, ml_path_node in zip(split_variants, ml_path_nodes):
                update_data = UpdateData(
                    ml_path_node_key=ml_path_node.key,
                    new_node_sequence=split_variant.get_mutated_sequence(ml_path_node)
                )
                update_data_list.append(update_data)

        return update_data_list


class DenovoVariantsDB:
    @staticmethod
    def _read_nb_of_samples(filehandler: TextIO) -> int:
        line = filehandler.readline().strip()
        nb_of_samples = int(line.split()[0])
        return nb_of_samples

    @staticmethod
    def _read_sample(filehandler: TextIO) -> str:
        line = filehandler.readline().strip()
        sample = line.split()[1]
        return sample

    @staticmethod
    def _read_nb_of_loci_in_sample(filehandler: TextIO) -> int:
        line = filehandler.readline().strip()
        nb_of_loci_in_sample = int(line.split()[0])
        return nb_of_loci_in_sample

    @staticmethod
    def _read_locus(filehandler: TextIO) -> str:
        locus = filehandler.readline().strip()
        return locus

    @staticmethod
    def _read_nb_of_nodes_in_ml_path(filehandler: TextIO) -> int:
        line = filehandler.readline().strip()
        nb_of_nodes_in_ml_path = int(line.split()[0])
        return nb_of_nodes_in_ml_path

    @classmethod
    def _read_MLPathNode(cls, filehandler: TextIO) -> MLPathNode:
        line = filehandler.readline().strip()
        try:
            matches = cls.ml_path_regex.search(line)
            start_index = int(matches.group(1))

            # TODO: fix pandora to give us non-inclusive end intervals instead
            end_index = int(matches.group(2)) + 1
            sequence = matches.group(3)
        except Exception as exc:
            assert False, f"Failed matching ML path regex to line: {line}\n" \
                          f"Exception: {str(exc)}"

        ml_path_node = MLPathNode(key=(start_index, end_index), sequence=sequence)
        return ml_path_node

    @classmethod
    def _read_ml_path(cls, filehandler: TextIO) -> MLPath:
        nb_of_nodes_in_ml_path = cls._read_nb_of_nodes_in_ml_path(filehandler)
        ml_path = []
        for _ in range(nb_of_nodes_in_ml_path):
            try:
                ml_path_node = cls._read_MLPathNode(filehandler)
                ml_path.append(ml_path_node)
            except EmptyMLPathSequence:
                # if the ML path node has empty sequence, it has empty interval, so we just ignore it,
                # as it is not amenable to updates
                pass

        return MLPath(ml_path)

    @staticmethod
    def _read_nb_of_variants(filehandler: TextIO) -> int:
        line = filehandler.readline().strip()
        nb_of_variants = int(line.split()[0])
        return nb_of_variants

    @staticmethod
    def _read_DenovoVariant(filehandler: TextIO) -> DenovoVariant:
        line = filehandler.readline()[:-1]
        line_split = line.split("\t")

        start_index_in_linear_path = int(line_split[0]) - 1
        ref = line_split[1]
        alt = line_split[2]
        denovo_variant = DenovoVariant(
                start_index_in_linear_path=start_index_in_linear_path,
                ref=ref,
                alt=alt,
        )
        return denovo_variant

    @classmethod
    def _read_variants(cls, filehandler) -> List[DenovoVariant]:
        nb_of_variants = cls._read_nb_of_variants(filehandler)
        variants = []
        for _ in range(nb_of_variants):
            denovo_variant = cls._read_DenovoVariant(filehandler)
            variants.append(denovo_variant)
        return variants

    def _get_locus_name_to_denovo_loci_core(self, filehandler: TextIO) -> Dict[str, List[DenovoLocusInfo]]:
        locus_name_to_denovo_loci = defaultdict(list)
        try:
            nb_of_samples = self._read_nb_of_samples(filehandler)
        except IndexError:
            logger.warning(
                f"File containing denovo paths ({self.filepath}) is empty, is it the correct file?"
            )
            return locus_name_to_denovo_loci

        for sample_index in range(nb_of_samples):
            sample = self._read_sample(filehandler)
            nb_of_loci_in_sample = self._read_nb_of_loci_in_sample(filehandler)
            for locus_index in range(nb_of_loci_in_sample):
                locus = self._read_locus(filehandler)
                ml_path = self._read_ml_path(filehandler)
                variants = self._read_variants(filehandler)
                denovo_locus = DenovoLocusInfo(sample, locus, ml_path, variants)
                locus_name_to_denovo_loci[locus].append(denovo_locus)

        return locus_name_to_denovo_loci

    def _get_locus_name_to_denovo_loci(self) -> Dict[str, List[DenovoLocusInfo]]:
        with open(self.filepath) as filehandler:
            return self._get_locus_name_to_denovo_loci_core(filehandler)

    @staticmethod
    def _get_locus_name_to_update_data(
            locus_name_to_denovo_loci: Dict[str, List[DenovoLocusInfo]])\
            -> Dict[str, List[UpdateData]]:
        locus_name_to_update_data = defaultdict(list)
        for locus_name, denovo_loci in locus_name_to_denovo_loci.items():
            for denovo_locus in denovo_loci:
                update_data = denovo_locus.get_update_data()
                locus_name_to_update_data[locus_name].extend(update_data)
        return locus_name_to_update_data

    def __init__(self, filepath: str):
        self.filepath: Path = Path(filepath)
        locus_name_to_denovo_loci = self._get_locus_name_to_denovo_loci()
        self.locus_name_to_update_data = self._get_locus_name_to_update_data(
            locus_name_to_denovo_loci
        )

    # Example:
    # (0 [0, 110) ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCACCGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC)
    ml_path_regex = re.compile(r"\(\d+ \[(\d+), (\d+)\) ([ACGT]*)\)")
