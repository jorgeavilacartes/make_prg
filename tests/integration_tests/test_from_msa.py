import os
from unittest import TestCase
from pathlib import Path
from make_prg.subcommands import from_msa
from argparse import Namespace
from tests.test_helpers import remove_dir_if_exists, are_dir_trees_equal
from make_prg.utils.seq_utils import SequenceCurationError


data_dir = Path("tests/integration_tests/data")


class Test_From_MSA_Integration_Full_Builds(TestCase):
    def prepare_options(self, test_name):
        input_data = str(data_dir / f"{test_name}.fa")
        output_folder = data_dir / "output" / test_name
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / test_name)

        options = Namespace(
            input=input_data,
            output_prefix=output_prefix,
            alignment_format='fasta',
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_graphs=False,
            threads=1,
            verbose=False)

        return options

    def test___match(self):
        options = self.prepare_options("match")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/match",
                                            data_dir / "output/match"))

    def test___match_nonmatch(self):
        options = self.prepare_options("match.nonmatch")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/match.nonmatch",
                                            data_dir / "output/match.nonmatch"))

    def test___match_nonmatch_match(self):
        options = self.prepare_options("match.nonmatch.match")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/match.nonmatch.match",
                                            data_dir / "output/match.nonmatch.match"))

    def test___match_nonmatch_shortmatch(self):
        options = self.prepare_options("match.nonmatch.shortmatch")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/match.nonmatch.shortmatch",
                                            data_dir / "output/match.nonmatch.shortmatch"))

    def test___match_staggereddash(self):
        options = self.prepare_options("match.staggereddash")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/match.staggereddash",
                                            data_dir / "output/match.staggereddash"))

    def test___nonmatch(self):
        options = self.prepare_options("nonmatch")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/nonmatch",
                                            data_dir / "output/nonmatch"))

    def test___nonmatch_match(self):
        options = self.prepare_options("nonmatch.match")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/nonmatch.match",
                                            data_dir / "output/nonmatch.match"))

    def test___nonmatch_shortmatch(self):
        options = self.prepare_options("nonmatch.shortmatch")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/nonmatch.shortmatch",
                                            data_dir / "output/nonmatch.shortmatch"))

    def test___shortmatch_nonmatch(self):
        options = self.prepare_options("shortmatch.nonmatch")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/shortmatch.nonmatch",
                                            data_dir / "output/shortmatch.nonmatch"))

    def test___shortmatch_nonmatch_match(self):
        options = self.prepare_options("shortmatch.nonmatch.match")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/shortmatch.nonmatch.match",
                                            data_dir / "output/shortmatch.nonmatch.match"))

    def test___contains_n(self):
        options = self.prepare_options("contains_n")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/contains_n",
                                            data_dir / "output/contains_n"))

    def test___contains_n_and_RYKMSW(self):
        options = self.prepare_options("contains_n_and_RYKMSW")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/contains_n_and_RYKMSW",
                                            data_dir / "output/contains_n_and_RYKMSW"))

    def test___contains_n_no_variants(self):
        options = self.prepare_options("contains_n_no_variants")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/contains_n_no_variants",
                                            data_dir / "output/contains_n_no_variants"))

    def test___contains_RYKMSW(self):
        options = self.prepare_options("contains_RYKMSW")
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/contains_RYKMSW",
                                            data_dir / "output/contains_RYKMSW"))

    def test___fails___a_column_full_of_Ns(self):
        options = self.prepare_options("fails")
        with self.assertRaises(SequenceCurationError):
            from_msa.run(options)

    def test___fails___unexpected_char_in_MSA(self):
        options = self.prepare_options("fails_2")
        with self.assertRaises(SequenceCurationError):
            from_msa.run(options)

    def test___nested_snp_backgrounds(self):
        options = self.prepare_options("nested_snps_seq_backgrounds")
        options.min_match_length = 3
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/nested_snps_seq_backgrounds",
                                            data_dir / "output/nested_snps_seq_backgrounds"))

    def test___nested_snps_under_del(self):
        options = self.prepare_options("nested_snps_deletion")
        options.min_match_length = 1
        from_msa.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/nested_snps_deletion",
                                            data_dir / "output/nested_snps_deletion"))

    def test___input_file_does_not_exist___raises_FileNotFoundError(self):
        options = self.prepare_options("unexistent_file")
        with self.assertRaises(FileNotFoundError):
            from_msa.run(options)

    # This tests the multi loci input and output
    def test___several_alignments(self):
        input_data = str(data_dir / "several")
        output_folder = data_dir / "output" / "several"
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / "several")

        options = Namespace(
            input=input_data,
            output_prefix=output_prefix,
            alignment_format='fasta',
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_graphs=False,
            threads=1,
            verbose=False)

        from_msa.run(options)

        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/several",
                                            data_dir / "output/several"))

    def test___input_dir_has_no_files___raises_FileNotFoundError(self):
        unexistent_folder = data_dir / "unexistent_folder"
        os.makedirs(unexistent_folder, exist_ok=True)
        input_data = str(unexistent_folder)
        output_folder = data_dir / "output" / "unexistent_folder"
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / "unexistent_folder")

        options = Namespace(
            input=input_data,
            output_prefix=output_prefix,
            alignment_format='fasta',
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_graphs=False,
            threads=1,
            verbose=False)

        with self.assertRaises(FileNotFoundError):
            from_msa.run(options)

    def test___output_files_already_exist(self):
        input_data = str(data_dir / "match.fa")
        output_prefix = str(data_dir / "truth_output/match/match")

        options = Namespace(
            input=input_data,
            output_prefix=output_prefix,
            alignment_format='fasta',
            log=None,
            max_nesting=5,
            min_match_length=7,
            output_graphs=False,
            threads=1,
            verbose=False)

        with self.assertRaises(RuntimeError):
            from_msa.run(options)