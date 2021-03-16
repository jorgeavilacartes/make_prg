from pathlib import Path
import logging
import sys
from datetime import datetime


def output_files_already_exist(output_prefix):
    return Path(output_prefix + "_prgs").exists() or \
           Path(output_prefix + "_tmp").exists() or \
           Path(output_prefix + ".prg.fa").exists() or \
           Path(output_prefix + ".update_DS").exists()


datefmt="%d/%m/%Y %I:%M:%S"
formatter = logging.Formatter(
    fmt="%(levelname)s %(asctime)s %(message)s", datefmt=datefmt
)


def setup_logging():
    handler_stderr = logging.StreamHandler(sys.stdout)
    handler_stderr.setFormatter(formatter)
    logging.getLogger().addHandler(handler_stderr)


def print_with_time(message):
    dateTimeObj = datetime.now()
    timestampStr = dateTimeObj.strftime(datefmt)
    print(f"{timestampStr} {message}")