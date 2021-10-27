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
        denovo_variant = DenovoVariant(5, "AGCT", "TGCA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "AG", "TG"), DenovoVariant(7, "CT", "CA")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___phased_snps_var_goes_through_two_nodes_split_in_start(self):
        denovo_variant = DenovoVariant(5, "ACGT", "TCTT")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2, dummy_node_2, dummy_node_2]

        expected = [DenovoVariant(5, "A", "T"), DenovoVariant(6, "CGT", "CTT")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___phased_snps_var_goes_through_two_nodes_split_in_end(self):
        denovo_variant = DenovoVariant(5, "ACGT", "AGGA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "ACG", "AGG"), DenovoVariant(8, "T", "A")]
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

    def test___split_variant___insertion_var_goes_through_two_nodes_split_in_start___one_variant_removed(self):
        denovo_variant = DenovoVariant(5, "TG", "CATG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "T", "CAT")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___insertion_var_goes_through_two_nodes_split_in_end___one_variant_removed(self):
        denovo_variant = DenovoVariant(5, "TG", "TGCA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(6, "G", "GCA")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)


    def test___split_variant___long_insertion_var_before_and_after_goes_through_two_nodes(self):
        denovo_variant = DenovoVariant(5, "TG", "AAAATGCCCC")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "T", "AAAAT"), DenovoVariant(6, "G", "GCCCC")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___long_insertion_var_between_goes_through_two_nodes(self):
        denovo_variant = DenovoVariant(5, "TG", "TAAAAG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(6, "G", "AAAAG")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___long_insertion_var_before_after_and_between_goes_through_two_nodes(self):
        denovo_variant = DenovoVariant(5, "TG", "AAAATCCCCGAAAA")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [dummy_node_1, dummy_node_2]

        expected = [DenovoVariant(5, "T", "AAAAT"), DenovoVariant(6, "G", "CCCCGAAAA")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___del_snp_ins_goes_through_three_nodes(self):
        denovo_variant = DenovoVariant(5, "AAAAAACGGGGGCTTTT", "AATGGGGGATTCCTTCC")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        dummy_node_3 = MLPathNode(key=(2, 3), sequence="T")
        # all 4 bases goes through this single node
        ml_path_nodes_it_goes_through = [
            dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_1, dummy_node_1,
            dummy_node_2, dummy_node_2, dummy_node_2, dummy_node_2, dummy_node_2, dummy_node_2, dummy_node_2,
            dummy_node_3, dummy_node_3, dummy_node_3, dummy_node_3]

        expected = [DenovoVariant(5, "AAAAAA", "AA"), DenovoVariant(11, "CGGGGGC", "TGGGGGA"),
                    DenovoVariant(18, "TTTT", "TTCCTTCC")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___single_SNPs_through_several_nodes(self):
        denovo_variant = DenovoVariant(5, "ACGTTCGT", "TCCAACGG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        dummy_node_3 = MLPathNode(key=(2, 3), sequence="T")
        dummy_node_4 = MLPathNode(key=(3, 4), sequence="A")
        dummy_node_5 = MLPathNode(key=(4, 5), sequence="C")
        dummy_node_6 = MLPathNode(key=(5, 6), sequence="T")
        dummy_node_7 = MLPathNode(key=(6, 7), sequence="A")
        dummy_node_8 = MLPathNode(key=(7, 8), sequence="C")
        ml_path_nodes_it_goes_through = [
            dummy_node_1, dummy_node_2, dummy_node_3, dummy_node_4, dummy_node_5, dummy_node_6, dummy_node_7,
            dummy_node_8]

        expected = [DenovoVariant(5, "A", "T"), DenovoVariant(7, "G", "C"), DenovoVariant(8, "T", "A"),
                    DenovoVariant(9, "T", "A"), DenovoVariant(12, "T", "G")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___split_variant___dels_through_several_nodes(self):
        denovo_variant = DenovoVariant(5, "ACGTTCGT", "AGTG")
        dummy_node_1 = MLPathNode(key=(0, 1), sequence="A")
        dummy_node_2 = MLPathNode(key=(1, 2), sequence="C")
        dummy_node_3 = MLPathNode(key=(2, 3), sequence="T")
        dummy_node_4 = MLPathNode(key=(3, 4), sequence="A")
        dummy_node_5 = MLPathNode(key=(4, 5), sequence="C")
        dummy_node_6 = MLPathNode(key=(5, 6), sequence="T")
        dummy_node_7 = MLPathNode(key=(6, 7), sequence="A")
        dummy_node_8 = MLPathNode(key=(7, 8), sequence="C")
        ml_path_nodes_it_goes_through = [
            dummy_node_1, dummy_node_2, dummy_node_3, dummy_node_4, dummy_node_5, dummy_node_6, dummy_node_7,
            dummy_node_8]

        expected = [DenovoVariant(6, "C", ""), DenovoVariant(8, "T", ""), DenovoVariant(10, "C", ""),
                    DenovoVariant(12, "T", "")]
        actual = denovo_variant.split_variant(ml_path_nodes_it_goes_through)

        self.assertEqual(expected, actual)

    def test___is_strict_insertion_event___true(self):
        denovo_variant = DenovoVariant(5, "", "AGTG")
        self.assertTrue(denovo_variant.is_strict_insertion_event())

    def test___is_strict_insertion_event___false(self):
        denovo_variant = DenovoVariant(5, "A", "AGTG")
        self.assertFalse(denovo_variant.is_strict_insertion_event())

    def test___str(self):
        expected = "[5:9]:'AACC'->'GT'"
        actual = str(self.sample_denovo_variant)
        self.assertEqual(expected, actual)

    def test___repr(self):
        expected = "DenovoVariant(start_index_in_linear_path=5, ref=\"AACC\", alt=\"GT\")"
        actual = repr(self.sample_denovo_variant)
        self.assertEqual(expected, actual)
