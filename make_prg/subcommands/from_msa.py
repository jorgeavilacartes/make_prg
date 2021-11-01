from typing import List
import multiprocessing
from pathlib import Path
from loguru import logger
from make_prg import prg_builder
from make_prg.from_msa import NESTING_LVL, MIN_MATCH_LEN
from make_prg.utils import io_utils, gfa


def register_parser(subparsers):
    subparser_msa = subparsers.add_parser(
        "from_msa",
        usage="make_prg from_msa",
        help="Make PRG from multiple sequence alignment",
    )
    subparser_msa.add_argument(
        "-i",
        "--input",
        action="store",
        type=str,
        required=True,
        help=(
            "Input data. If it is a file, will be read as the supported alignment_format (see "
            "parameter --alignment_format). A single PRG, GFA and update data structures will "
            "be created for the given MSA. If it is a dir, will process all files in the dir. "
            "Multiple PRGs, GFAs and update data structures will be created, one for each MSA"
        ),
    )
    subparser_msa.add_argument(
        "-o",
        "--output_prefix",
        action="store",
        type=str,
        required=True,
        help="Prefix for the output files",
    )
    subparser_msa.add_argument(
        "-t",
        "--threads",
        action="store",
        type=int,
        default=1,
        help="Number of threads (default 1)",
    )
    subparser_msa.add_argument(
        "-f",
        "--alignment_format",
        dest="alignment_format",
        action="store",
        default="fasta",
        help=(
            "Alignment format of MSA, must be a biopython AlignIO input "
            "alignment_format. See http://biopython.org/wiki/AlignIO. Default: %(default)s"
        ),
    )
    subparser_msa.add_argument(
        "-N",
        "--max_nesting",
        dest="max_nesting",
        action="store",
        type=int,
        default=NESTING_LVL,
        help="Maximum number of levels to use for nesting. Default: %(default)s",
    )
    subparser_msa.add_argument(
        "-L",
        "--min_match_length",
        dest="min_match_length",
        action="store",
        type=int,
        default=MIN_MATCH_LEN,
        help=(
            "Minimum number of consecutive characters which must be identical for a "
            "match. Default: %(default)s"
        ),
    )
    subparser_msa.add_argument(
        "--output_graphs",
        dest="output_graphs",
        action="store_true",
        default=False,
        help="Outputs the recursive tree graphical representation (for development use only)",
    )

    subparser_msa.set_defaults(func=run)

    return subparser_msa


def get_all_input_files(input_path: str) -> List[Path]:
    input_path = Path(input_path)
    if not input_path.exists():
        raise FileNotFoundError(f"{input_path} does not exist")

    if input_path.is_file():
        all_files = [input_path]
    else:
        all_files = [
            Path(path).absolute() for path in input_path.iterdir() if path.is_file()
        ]
    return all_files


def process_MSA(msa_filepath: Path):
    logger.info(f"Generating PRG for {msa_filepath}...")
    msa_name = msa_filepath.name
    locus_name = msa_filepath.with_suffix("").name

    temp_dir = io_utils.get_temp_dir_for_multiprocess(Path(options.output_prefix))
    prefix = str(temp_dir / msa_name)

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
        builder.write_prg(prefix, prg)
        gfa.GFA_Output.write_gfa(prefix, prg)
        builder.serialize(f"{prefix}.pickle")
        if options.output_graphs:
            builder.output_debug_graphs(Path(options.output_prefix + "_debug_graphs"))

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
        raise FileNotFoundError(f"No input files found in {options.input}")

    is_a_single_file = len(input_files) == 1
    if is_a_single_file:
        logger.info("A single file was given as input, single outputs will be produced")
    else:
        logger.info("Multiple files were given as input, multiple outputs will be produced")

    if io_utils.output_files_already_exist(options.output_prefix):
        raise RuntimeError("One or more output files already exists, aborting run...")

    logger.info(f"Using {options.threads} threads to generate PRGs...")
    with multiprocessing.Pool(options.threads) as pool:
        pool.map(process_MSA, input_files, chunksize=1)
    logger.success(f"All PRGs generated!")

    io_utils.create_final_files(options.threads, options.output_prefix,
                                is_a_single_MSA=is_a_single_file)
    logger.success("All done!")
