from unittest import TestCase
from pathlib import Path
from make_prg.subcommands import update
from argparse import Namespace
from tests.test_helpers import remove_dir_if_exists, are_dir_trees_equal


data_dir = Path("tests/integration_tests/data")


# NOTE: for these tests we need mafft in $PATH
class Test_Update_Integration_Full_Builds(TestCase):
    def prepare_options(self, test_name: str, update_DS: Path):
        output_folder = data_dir / "output_update" / test_name
        remove_dir_if_exists(output_folder)
        output_prefix = str(output_folder / test_name)

        options = Namespace(
            denovo_paths=str(data_dir / test_name / 'denovo_paths.txt'),
            update_DS=update_DS,
            output_prefix=output_prefix,
            mafft='mafft',
            log=None,
            output_graphs=False,
            threads=1,
            verbose=False)

        return options

    def test___update___match_update_simple(self):
        options = self.prepare_options(test_name="match_update_simple",
                                       update_DS=data_dir/"truth_output/match/match.update_DS.zip")
        update.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output_update/match_update_simple",
                                            data_dir / "output_update/match_update_simple"))

    def test___update___match_update_complex___several_samples_and_variants(self):
        options = self.prepare_options(test_name="match_update_complex",
                                       update_DS=data_dir/"truth_output/match/match.update_DS.zip")
        update.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output_update/match_update_complex",
                                            data_dir / "output_update/match_update_complex"))

    def test___update___match_nonmatch_match___locus_with_different_ML_path_for_each_sample(self):
        options = self.prepare_options(test_name="match.nonmatch.match_update",
                                       update_DS=data_dir/"truth_output/match.nonmatch.match/match.nonmatch.match.update_DS.zip")
        update.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output_update/match.nonmatch.match_update",
                                            data_dir / "output_update/match.nonmatch.match_update"))

    # TODO: understand why sometimes fails and most of times is ok
    def test___update___nested_snp_backgrounds(self):
        options = self.prepare_options(test_name="nested_snps_seq_backgrounds_update",
                                       update_DS=data_dir/"truth_output/nested_snps_seq_backgrounds/nested_snps_seq_backgrounds.update_DS.zip")
        update.run(options)
        self.assertTrue(are_dir_trees_equal(data_dir / "truth_output_update/nested_snps_seq_backgrounds_update",
                                            data_dir / "output_update/nested_snps_seq_backgrounds_update"))

    # TODO: finish
    # def test___update___nested_snps_under_del(self):
    #     options = self.prepare_options("nested_snps_deletion")
    #     options.min_match_length = 1
    #     update.run(options)
    #     self.assertTrue(are_dir_trees_equal(data_dir / "truth_output/nested_snps_deletion",
    #                                         data_dir / "output/nested_snps_deletion"))

    def test___update___output_files_already_exist(self):
        output_folder = data_dir / "output_update" / "match_update_simple"
        output_prefix = str(output_folder / "match_update_simple")

        options = Namespace(
            denovo_paths=str(data_dir / "match_update_simple" / 'denovo_paths.txt'),
            update_DS=data_dir / "truth_output/match/match.update_DS.zip",
            output_prefix=output_prefix,
            mafft='mafft',
            log=None,
            output_graphs=False,
            threads=1,
            verbose=False)

        with self.assertRaises(RuntimeError):
            update.run(options)
