import logging
import os
from pathlib import Path

from make_prg import io_utils
from make_prg.prg_builder import PrgBuilder
from make_prg.utils import output_files_already_exist
from make_prg.denovo_paths_reader import DenovoPathsDB

def register_parser(subparsers):
    subparser_update_prg = subparsers.add_parser(
        "update",
        usage="make_prg update_prg",
        help="Update PRGs given new sequences output by pandora.",
    )
    subparser_update_prg.add_argument(
        "-p", "--pickle",
        action="store",
        type=str,
        required=True,
        help=(
            "Pickle filepath containing PRGs created previously with this tool."
        ),
    )
    subparser_update_prg.add_argument(
        "-d", "--denovo_paths",
        action="store",
        type=str,
        required=True,
        help=(
            "Filepath containing denovo sequences output by pandora (denovo_paths.txt)."
        ),
    )
    subparser_update_prg.add_argument(
        "-o", "--output_prefix",
        action="store",
        type=str,
        required=True,
        help=(
            "Output prefix: prefix for the output files"
        ),
    )
    subparser_update_prg.add_argument(
        "-t", "--threads",
        action="store",
        type=int,
        default=1,
        help="Number of threads",
    )
    subparser_update_prg.set_defaults(func=run)

    return subparser_update_prg


def update(locus_name_to_prg_builder, denovo_paths_db, output_prefix, threads):
    # dummy entrypoint for now
    print("Update feature called")


def run(options):
    if output_files_already_exist(options.output_prefix):
        raise RuntimeError("One or more output files already exists, aborting run...")

    locus_name_to_prg_builder = PrgBuilder.deserialize(options.pickle)
    denovo_paths_db = DenovoPathsDB(options.denovo_paths)
    output_dir = Path(options.output_prefix).parent
    os.makedirs(output_dir, exist_ok=True)
    update(locus_name_to_prg_builder, denovo_paths_db, options.output_prefix, options.threads)