from typing import List
import multiprocessing
import os
from pathlib import Path
from loguru import logger
from make_prg import prg_builder
from make_prg.from_msa import NESTING_LVL, MIN_MATCH_LEN
from make_prg.utils import io_utils
from make_prg.utils.utils import output_files_already_exist
from make_prg.utils.drawer import RecursiveTreeDrawer


def register_parser(subparsers):
    subparser_msa = subparsers.add_parser(
        "from_msa",
        usage="make_prg from_msa",
        help="Make PRG from multiple sequence alignment dir",
    )
    subparser_msa.add_argument(
        "-i",
        "--input",
        action="store",
        type=str,
        required=True,
        help=(
            "Input dir: all files in this will try to be read as the supported alignment_format. "
            "If not aligned in fasta alignment_format, use -f to input the "
            "alignment_format type"
        ),
    )
    subparser_msa.add_argument(
        "-o",
        "--output_prefix",
        action="store",
        type=str,
        required=True,
        help=("Output prefix: prefix for the output files"),
    )
    subparser_msa.add_argument(
        "-t",
        "--threads",
        action="store",
        type=int,
        default=1,
        help="Number of threads",
    )
    subparser_msa.add_argument(
        "-f",
        "--alignment_format",
        dest="alignment_format",
        action="store",
        default="fasta",
        help=(
            "Alignment format of MSA, must be a biopython AlignIO input "
            "alignment_format. See http://biopython.org/wiki/AlignIO. Default: fasta"
        ),
    )
    subparser_msa.add_argument(
        "--max_nesting",
        dest="max_nesting",
        action="store",
        type=int,
        default=NESTING_LVL,
        help="Maximum number of levels to use for nesting. Default: {}".format(
            NESTING_LVL
        ),
    )
    subparser_msa.add_argument(
        "--min_match_length",
        dest="min_match_length",
        action="store",
        type=int,
        default=MIN_MATCH_LEN,
        help=(
            "Minimum number of consecutive characters which must be identical for a "
            "match. Default: {}".format(MIN_MATCH_LEN)
        ),
    )
    subparser_msa.add_argument(
        "--output_graphs",
        dest="output_graphs",
        action="store_true",
        default=False,
        help="Outputs the recursive tree and the PRG graphical representation (for development use only)",
    )
    subparser_msa.set_defaults(func=run)

    return subparser_msa


def get_all_input_files(input_dir: str) -> List[Path]:
    input_dir = Path(input_dir)
    all_files = [
        Path(path).absolute() for path in input_dir.iterdir() if path.is_file()
    ]
    return all_files


def process_MSA(msa_filepath: Path):
    logger.info(f"Generating PRG for {msa_filepath}...")
    msa_name = msa_filepath.name
    locus_name = msa_filepath.with_suffix("").name
    current_process = multiprocessing.current_process()

    workdir = Path(options.output_prefix + "_prgs") / current_process.name
    os.makedirs(workdir, exist_ok=True)
    prefix = str(workdir / msa_name)

    if options.output_graphs:
        debug_graphs_dir = Path(options.output_prefix + "_debug_graphs")
        os.makedirs(debug_graphs_dir, exist_ok=True)

    try:
        builder = prg_builder.PrgBuilder(
            locus_name=locus_name,
            msa_file=msa_filepath,
            alignment_format=options.alignment_format,
            max_nesting=options.max_nesting,
            min_match_length=options.min_match_length
        )

        logger.info(f"Writing output files of locus {msa_name}")
        prg = builder.build_prg()
        io_utils.write_prg(prefix, prg)
        builder.serialize(f"{prefix}.pickle")

        if options.output_graphs:
            recursive_tree_drawer = RecursiveTreeDrawer(builder.root)
            recursive_tree_drawer.output_graph(debug_graphs_dir / f"{msa_name}.recursion_tree.png")

        # TODO: add back GFA writing
        # print_with_time(f"Write GFA file to {prefix}.gfa")
        # io_utils.write_gfa(f"{prefix}.gfa", aseq.prg)
    except ValueError as value_error:
        if "No records found in handle" in value_error.args[0]:
            logger.warning(f"No records found in MSA {msa_filepath}, skipping...")
        else:
            raise value_error


def run(cl_options):
    global options
    options = cl_options

    logger.info("Getting input files...")
    input_files = get_all_input_files(options.input)
    there_is_no_input_files = len(input_files) == 0
    if there_is_no_input_files:
        logger.warning(f"No input files found at {options.input}")

    if output_files_already_exist(options.output_prefix):
        raise RuntimeError("One or more output files already exists, aborting run...")

    # NB: don't use logging, it causes deadlocks: https://pythonspeed.com/articles/python-multiprocessing/
    logger.info(f"Using {options.threads} threads to generate PRGs...")
    with multiprocessing.Pool(options.threads) as pool:
        pool.map(process_MSA, input_files, chunksize=1)
    logger.success(f"All PRGs generated!")

    logger.info("Concatenating files from several threads into single final files...")
    # get all files that were generated
    prg_files = []
    locus_name_to_pickle_files = {}
    for process_num in range(1, options.threads + 1):
        workdir = (
            Path(options.output_prefix + "_prgs") / f"ForkPoolWorker-{process_num}"
        )
        if workdir.exists():
            for file in workdir.iterdir():
                if file.is_file():
                    if file.name.endswith(".prg.fa"):
                        prg_files.append(file)
                    elif file.name.endswith(".pickle"):
                        locus_name = file.with_suffix("").with_suffix("").name
                        locus_name_to_pickle_files[locus_name] = file

    # concatenate the prg.fa output files
    io_utils.concatenate_text_files(prg_files, options.output_prefix + ".prg.fa")
    # cleanup
    for prg_file in prg_files:
        prg_file.unlink()

    # create and serialise the PRG Builder collection
    logger.info("Serialising update data structure...")
    prg_builder_collection = prg_builder.PrgBuilderCollection(Path(f"{options.output_prefix}.update_DS.zip"))
    prg_builder_collection.save(locus_name_to_pickle_files)

    # cleanup
    for pickle_file in locus_name_to_pickle_files.values():
        pickle_file.unlink()
    io_utils.remove_empty_folders(options.output_prefix + "_prgs")


    logger.success("All done!")
