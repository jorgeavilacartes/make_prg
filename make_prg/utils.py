from pathlib import Path


def output_files_already_exist(output_prefix):
    return (
        Path(output_prefix + "_prgs").exists()
        or Path(output_prefix + "_tmp").exists()
        or Path(output_prefix + ".prg.fa").exists()
        or Path(output_prefix + ".update_DS").exists()
    )


def assert_sequence_is_composed_of_ACGT_only(seq):
    for base in seq:
        assert base in "ACGT"
