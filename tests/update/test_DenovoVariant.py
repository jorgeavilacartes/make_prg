from unittest import TestCase
from make_prg.update.denovo_variants import DenovoVariant, DenovoError
from make_prg.update.MLPath import MLPathNode

class DenovoVariantTest(TestCase):
    def setUp(self) -> None:
        self.sample_denovo_variant = DenovoVariant(5, "AACC", "GT")

    def test___ref_not_composed_of_ACGT_only___DenovoError_raised(self):
        with self.assertRaises(DenovoError):
            DenovoVariant(0, ref="ACGTB", alt="")

    def test___alt_not_composed_of_ACGT_only___DenovoError_raised(self):
        with self.assertRaises(DenovoError):
            DenovoVariant(0, ref="", alt="ACGTB")

    def test___ref_and_alt_are_identical___DenovoError_raised(self):
        seq="ACGT"
        with self.assertRaises(DenovoError):
            DenovoVariant(0, ref=seq, alt=seq)

    def test___negative_index_for_variant_pos___DenovoError_raised(self):
        with self.assertRaises(DenovoError):
            DenovoVariant(-1, "A", "C")

    def test___constructor___variant_correctly_built(self):
        self.assertEqual(5, self.sample_denovo_variant.start_index_in_linear_path)
        self.assertEqual(9, self.sample_denovo_variant.end_index_in_linear_path)
        self.assertEqual("AACC", self.sample_denovo_variant.ref)
        self.assertEqual("GT", self.sample_denovo_variant.alt)

    def test___get_mutated_sequence___bad_start_index_before_node_start_DenovoError_raised(self):
        ml_path_node = MLPathNode(key=(6, 9), sequence="ACG")
        ml_path_node.start_index_in_linear_path = 6
        ml_path_node.end_index_in_linear_path = 9

        with self.assertRaises(DenovoError):
            self.sample_denovo_variant.get_mutated_sequence(ml_path_node)

    def test___get_mutated_sequence___bad_end_index_after_node_end_DenovoError_raised(self):
        ml_path_node = MLPathNode(key=(5, 8), sequence="ACG")
        ml_path_node.start_index_in_linear_path = 5
        ml_path_node.end_index_in_linear_path = 8

        with self.assertRaises(DenovoError):
            self.sample_denovo_variant.get_mutated_sequence(ml_path_node)

    def test___get_mutated_sequence___node_seq_not_consistent_with_ref_DenovoError_raised(self):
        ml_path_node = MLPathNode(key=(3, 13), sequence="GTAAGCGTCA")
        ml_path_node.start_index_in_linear_path = 3
        ml_path_node.end_index_in_linear_path = 13

        with self.assertRaises(DenovoError):
            self.sample_denovo_variant.get_mutated_sequence(ml_path_node)

    def test___get_mutated_sequence(self):
        ml_path_node = MLPathNode(key=(3, 13), sequence="GTAACCGTCA")
        ml_path_node.start_index_in_linear_path = 3
        ml_path_node.end_index_in_linear_path = 13

        expected = "GTGTGTCA"
        actual = self.sample_denovo_variant.get_mutated_sequence(ml_path_node)

        self.assertEqual(expected, actual)

    def test___split_variant___empty_path_error_out(self):
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.split_variant([])

    def test___split_variant___one_node_less_than_expected_in_path_error_out(self):
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        ml_path_nodes_it_goes_through = [dummy_node, dummy_node, dummy_node]
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.split_variant(ml_path_nodes_it_goes_through)

    def test___split_variant___one_node_more_than_expected_in_path_error_out(self):
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        ml_path_nodes_it_goes_through = [dummy_node, dummy_node, dummy_node, dummy_node, dummy_node]
        with self.assertRaises(AssertionError):
            self.sample_denovo_variant.split_variant(ml_path_nodes_it_goes_through)

    def test___split_variant___variant_goes_through_single_node___no_split(self):
        dummy_node = MLPathNode(key=(0, 1), sequence="A")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node, dummy_node, dummy_node, dummy_node]

        expected = [self.sample_denovo_variant]
        actual = self.sample_denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_middle(self):
        denovo_variant = DenovoVariant(5, "ACGT", "CG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AC", "C"), DenovoVariant(7, "GT", "G")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_variant_goes_through_two_nodes_split_in_middle_with_mismatch(self):
        denovo_variant = DenovoVariant(5, "ACGTA", "CAT")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "ACG", "CA"), DenovoVariant(8, "TA", "T")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_start(self):
        denovo_variant = DenovoVariant(5, "ACGT", "GCG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AC", "GC"), DenovoVariant(7, "GT", "G")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_start___one_variant_removed(self):
        denovo_variant = DenovoVariant(5, "ACGT", "AC")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(7, "GT", "")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_end(self):
        denovo_variant = DenovoVariant(5, "ACGT", "CCT")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AC", "C"), DenovoVariant(7, "GT", "CT")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___deletion_var_goes_through_two_nodes_split_in_end___one_variant_removed(self):
        denovo_variant = DenovoVariant(5, "ACGT", "GT")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AC", "")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___phased_snps_var_goes_through_two_nodes_split_in_middle(self):
        denovo_variant = DenovoVariant(5, "ACGT", "TGCA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AC", "TG"), DenovoVariant(7, "GT", "CA")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___phased_snps_var_goes_through_two_nodes_split_in_start(self):
        denovo_variant = DenovoVariant(5, "ACGT", "TGCA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "A", "T"), DenovoVariant(6, "CGT", "GCA")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___phased_snps_var_goes_through_two_nodes_split_in_end(self):
        denovo_variant = DenovoVariant(5, "ACGT", "TGCA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "ACG", "TGC"), DenovoVariant(8, "T", "A")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___insertion_var_goes_through_two_nodes_split_in_middle(self):
        denovo_variant = DenovoVariant(5, "GC", "TGCA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "G", "TG"), DenovoVariant(6, "C", "CA")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

