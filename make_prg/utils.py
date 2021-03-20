from pathlib import Path
import logging
import sys
from datetime import datetime


def output_files_already_exist(output_prefix):
    return Path(output_prefix + "_prgs").exists() or \
           Path(output_prefix + "_tmp").exists() or \
           Path(output_prefix + ".prg.fa").exists() or \
           Path(output_prefix + ".update_DS").exists()


datefmt="%d/%m/%Y %H:%M:%S"
formatter = logging.Formatter(
    fmt="%(levelname)s %(asctime)s %(message)s", datefmt=datefmt
)


def setup_logging():
    handler_stderr = logging.StreamHandler(sys.stdout)
    handler_stderr.setFormatter(formatter)
    logging.getLogger().addHandler(handler_stderr)


def print_with_time(message):
    date_time_obj = datetime.now()
    timestamp_str = date_time_obj.strftime(datefmt)
    print(f"{timestamp_str} {message}")


def assert_sequence_is_composed_of_ACGT_only(seq):
    for base in seq:
        assert base in "ACGT"
