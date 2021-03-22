from pathlib import Path
import sys
from datetime import datetime


def output_files_already_exist(output_prefix):
    return Path(output_prefix + "_prgs").exists() or \
           Path(output_prefix + "_tmp").exists() or \
           Path(output_prefix + ".prg.fa").exists() or \
           Path(output_prefix + ".update_DS").exists()


datefmt="%d/%m/%Y %H:%M:%S"
def print_with_time(message):
    date_time_obj = datetime.now()
    timestamp_str = date_time_obj.strftime(datefmt)
    print(f"{timestamp_str} {message}")


def assert_sequence_is_composed_of_ACGT_only(seq):
    for base in seq:
        assert base in "ACGT"
