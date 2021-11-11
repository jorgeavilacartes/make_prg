from unittest import TestCase
from pathlib import Path
from make_prg.subcommands import from_msa
from argparse import Namespace
from tests.test_helpers import remove_dir_if_exists, are_dir_trees_equal


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

    # def test_answers_non_nested(self):
    #     infile = data_dir / "nonmatch.fa"
    #     aseq = PrgBuilder(infile)
    #     self.assertEqual(aseq.prg, " 5 AAACGTGGTT 6 CCCCCCCCCC 5 ")
    #
    #     infile = data_dir / "nonmatch.match.fa"
    #     aseq = PrgBuilder(infile)
    #     self.assertEqual(aseq.prg, " 5 AAACGT 6 CCCCCC 5 GGTT")
    #
    #     infile = data_dir / "shortmatch.nonmatch.match.fa"
    #     aseq = PrgBuilder(infile)
    #     self.assertEqual(aseq.prg, " 5 AAACGT 6 ATTTTC 5 GGTT")
    #
    #
    #     infile = data_dir / "contains_n.fa"
    #     aseq = PrgBuilder(infile)
    #     self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")
    #
    #     infile = data_dir / "contains_RYKMSW.fa"
    #     aseq = PrgBuilder(infile)
    #     self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")
    #
    #     infile = data_dir / "contains_n_and_RYKMSW.fa"
    #     aseq = PrgBuilder(infile)
    #     self.assertEqual(aseq.prg, "AAACG 5 T 6 C 5 GGTT")
    #
    #     infile = data_dir / "contains_n_and_RYKMSW_no_variants.fa"
    #     aseq = PrgBuilder(infile)
    #     self.assertEqual(aseq.prg, "AAACGTGGTT")
    #
    #     # with pytest.raises(Exception):
    #     #    aseq = AlignedSeq("test/fails.fa")
    #
    # def test_nested_snp_backgrounds(self):
    #     infile = data_dir / "nested_snps_seq_backgrounds.fa"
    #     aseq = PrgBuilder(infile, min_match_length=3)
    #     self.assertEqual(
    #         aseq.prg, " 5 AAAA 7 T 8 C 7 AAAAAA 6 CCCC 9 T 10 G 9 CCCCCC 5 "
    #     )
    #
    # def test_nested_snps_under_del(self):
    #     infile = data_dir / "nested_snps_deletion.fa"
    #     aseq = PrgBuilder(infile, min_match_length=1)
    #     self.assertEqual(aseq.prg, "A 5 AA 7 C 8 T 7 AAAA 9 T 10 G 9 AA 6 A 5 AA")
