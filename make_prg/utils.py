from pathlib import Path
import logging
import os
import sys

def output_files_already_exist(output_prefix):
    return Path(output_prefix + "_prgs").exists() or \
           Path(output_prefix + "_tmp").exists() or \
           Path(output_prefix + ".prg.fa").exists() or \
           Path(output_prefix + ".update_DS").exists()


formatter = logging.Formatter(
    fmt="%(levelname)s %(asctime)s %(message)s", datefmt="%d/%m/%Y %I:%M:%S"
)

# TODO: this might not work for multiprocesses...
def setup_file_logging(prefix):
    log_file = f"{prefix}.log"
    if os.path.exists(log_file):
        os.unlink(log_file)
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)
    return handler

def setup_stderr_logging():
    handler_stderr = logging.StreamHandler(sys.stderr)
    handler_stderr.setFormatter(formatter)
    logging.getLogger().addHandler(handler_stderr)