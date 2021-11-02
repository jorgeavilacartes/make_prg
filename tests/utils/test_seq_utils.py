from unittest import TestCase

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from tests.from_msa import make_alignment
from make_prg.utils.seq_utils import (
    is_non_match,
    NONMATCH,
    is_gap,
    GAP,
    remove_duplicates,
    ungap,
    get_number_of_unique_ungapped_sequences,
    get_number_of_unique_gapped_sequences,
    get_alignment_seqs,
    count,
    get_expanded_sequences,
    has_empty_sequence,
)


class TestSeqUtilsMisc(TestCase):
    def test___is_non_match___column_is_all_As(self):
        expected = False
        actual = is_non_match("A")
        self.assertEqual(expected, actual)

    def test___is_non_match___column_is_non_match(self):
        expected = True
        actual = is_non_match(NONMATCH)
        self.assertEqual(expected, actual)

    def test___is_gap___column_is_all_As(self):
        expected = False
        actual = is_gap("A")
        self.assertEqual(expected, actual)

    def test___is_gap___column_is_gap(self):
        expected = True
        actual = is_gap(GAP)
        self.assertEqual(expected, actual)

    def test___ungap___seq_has_no_gaps(self):
        sequence = "ACGTGTGACA"
        expected = "ACGTGTGACA"
        actual = ungap(sequence)
        self.assertEqual(expected, actual)

    def test___ungap___seq_has_gaps(self):
        sequence = "---AC-GT---GTG--AC--A"
        expected = "ACGTGTGACA"
        actual = ungap(sequence)
        self.assertEqual(expected, actual)

    def test___get_number_of_unique_ungapped_sequences___single_spaces(self):
        msa = make_alignment([
            "-ACGT",
            "A-CGT",
            "AC-GT",
            "ACG-T",
            "ACGT-",
            "AC--T",
            "A-CGT",
            "A-CGT",
            "ACG-T",
            "ACG-T",
            "AC--T",
            "AC--T",
            "AC--T",
            "A-CGT",
        ])

        expected = 2
        actual = get_number_of_unique_ungapped_sequences(msa)

        self.assertEqual(expected, actual)

    def test___get_number_of_unique_ungapped_sequences___several_spaces(self):
        msa = make_alignment([
            "TA--------T",
            "T-A-------T",
            "T--A------T",
            "T---A-----T",
            "T----A----T",
            "T-----A---T",
            "T------A--T",
            "T-------A-T",
            "T--------AT",
            "-T-------AT",
            "--T------AT",
            "---T-----AT",
            "T-------AT-",
            "T------AT--",
            "T-----AT---",
            "T--A------T",
            "T--A------T",
            "---T-----AT",
            "---T-----AT",
            "--T------AT",
            "-T-------AT",
            "T-----A---T",
            "T------A--T",
            "T-------A-T",
        ])

        expected = 1
        actual = get_number_of_unique_ungapped_sequences(msa)

        self.assertEqual(expected, actual)

    def test___get_number_of_unique_gapped_sequences___single_spaces(self):
        msa = make_alignment([
            "-ACGT",
            "A-CGT",
            "AC-GT",
            "ACG-T",
            "ACGT-",
            "AC--T",
            "A-CGT",
            "A-CGT",
            "ACG-T",
            "ACG-T",
            "AC--T",
            "AC--T",
            "AC--T",
            "A-CGT",
        ])

        expected = 6
        actual = get_number_of_unique_gapped_sequences(msa)

        self.assertEqual(expected, actual)

    def test___get_number_of_unique_gapped_sequences___several_spaces(self):
        msa = make_alignment([
            "TA--------T",
            "T-A-------T",
            "T--A------T",
            "T---A-----T",
            "T----A----T",
            "T-----A---T",
            "T------A--T",
            "T-------A-T",
            "T--------AT",
            "-T-------AT",
            "--T------AT",
            "---T-----AT",
            "T-------AT-",
            "T------AT--",
            "T-----AT---",
            "T--A------T",
            "T--A------T",
            "---T-----AT",
            "---T-----AT",
            "--T------AT",
            "-T-------AT",
            "T-----A---T",
            "T------A--T",
            "T-------A-T",
        ])

        expected = 15
        actual = get_number_of_unique_gapped_sequences(msa)

        self.assertEqual(expected, actual)

class TestRemoveDuplicates(TestCase):
    def test___remove_duplicates___empty_list(self):
        the_list = []

        expected = []
        actual = list(remove_duplicates(the_list))

        self.assertEqual(expected, actual)

    def test___remove_duplicates___unique_list(self):
        the_list = ["1", "2", "3"]

        expected = ["1", "2", "3"]
        actual = list(remove_duplicates(the_list))

        self.assertEqual(expected, actual)

    def test___remove_duplicates___several_copies(self):
        the_list = ["1", "1", "2", "3", "2", "2", "2", "2", "1", "1"]

        expected = ["1", "2", "3"]
        actual = list(remove_duplicates(the_list))

        self.assertEqual(expected, actual)


class TestEmptySeq(TestCase):
    def test___has_empty_sequence___individual_columns(self):
        msa = make_alignment(["ACCCT", "A---T", "A-G-T"])
        self.assertFalse(has_empty_sequence(msa, (0, 0)))
        self.assertTrue(has_empty_sequence(msa, (1, 1)))
        self.assertTrue(has_empty_sequence(msa, (2, 2)))
        self.assertTrue(has_empty_sequence(msa, (3, 3)))
        self.assertFalse(has_empty_sequence(msa, (4, 4)))

    def test___has_empty_sequence___long_columns(self):
        msa = make_alignment(["AAACCCGGGTTT", "AAA-C----T--", "AAA----G---T"])
        self.assertFalse(has_empty_sequence(msa, (0, 2)))
        self.assertTrue(has_empty_sequence(msa, (3, 5)))
        self.assertTrue(has_empty_sequence(msa, (6, 8)))
        self.assertFalse(has_empty_sequence(msa, (9, 11)))


class TestGetExpandedSequences(TestCase):
    def test_ambiguous_bases_one_seq(self):
        alignment = AlignIO.MultipleSeqAlignment([SeqRecord(Seq("RWAAT"))])
        result = get_expanded_sequences(alignment)
        expected = {"GAAAT", "AAAAT", "GTAAT", "ATAAT"}
        self.assertEqual(set(result), expected)

    def test_ambiguous_bases_one_seq_with_repeated_base(self):
        alignment = AlignIO.MultipleSeqAlignment([SeqRecord(Seq("RRAAT"))])
        result = get_expanded_sequences(alignment)
        expected = {"GAAAT", "AAAAT", "GGAAT", "AGAAT"}
        self.assertEqual(set(result), expected)

    def test_first_sequence_in_is_first_sequence_out(self):
        alignment = make_alignment(["TTTT", "AAAA", "CC-C"])
        result = get_expanded_sequences(alignment)
        expected = ["TTTT", "AAAA", "CCC"]
        self.assertEqual(expected, result)


class TestSeqIteration(TestCase):
    input_seqs = ["TTAGGTTT", "TTA--TTT", "GGA-TTTT"]

    def test_get_seqs_from_alignment(self):
        msa = make_alignment(self.input_seqs)
        result = list(get_alignment_seqs(msa))
        self.assertEqual(result, self.input_seqs)

    def test_count_alignment_seqs(self):
        msa = make_alignment(self.input_seqs)
        result = count(get_alignment_seqs(msa))
        self.assertEqual(result, 3)