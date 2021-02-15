import logging
import os
from pathlib import Path

from make_prg.from_msa import prg_builder, NESTING_LVL, MIN_MATCH_LEN
from make_prg import io_utils
from multiprocessing import Pool

options = None

def register_parser(subparsers):
    subparser_msa = subparsers.add_parser(
        "from_msa",
        usage="make_prg from_msa [options] <MSA input dir>",
        help="Make PRG from multiple sequence alignment dir",
    )
    subparser_msa.add_argument(
        "-i", "--input",
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
        "-o", "--output",
        action="store",
        type=str,
        required=True,
        help=(
            "Output dir: where all output will be stored"
        ),
    )
    subparser_msa.add_argument(
        "-t", "--threads",
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
    subparser_msa.set_defaults(func=run)

    return subparser_msa


def get_all_input_files(input_dir):
    input_dir = Path(input_dir)
    all_files = [Path(path).absolute() for path in input_dir.iterdir() if path.is_file()]
    return all_files


def process_MSA(msa_filepath: Path):
    msa_name = msa_filepath.name
    output = Path(options.output)
    prefix = str(output / msa_name)
    prefix += ".max_nest%d.min_match%d" % (
        options.max_nesting,
        options.min_match_length,
    )

    # Set up file logging
    log_file = f"{prefix}.log"
    if os.path.exists(log_file):
        os.unlink(log_file)
    formatter = logging.Formatter(
        fmt="%(levelname)s %(asctime)s %(message)s", datefmt="%d/%m/%Y %I:%M:%S"
    )
    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)
    logging.getLogger().addHandler(handler)

    logging.info(
        "Input parameters max_nesting: %d, min_match_length: %d",
        options.max_nesting,
        options.min_match_length,
    )

    builder = prg_builder.PrgBuilder(
        locus_name=msa_name,
        msa_file=msa_filepath,
        alignment_format=options.alignment_format,
        max_nesting=options.max_nesting,
        min_match_length=options.min_match_length,
    )
    builder.build_prg()
    logging.info(f"Write PRG file to {prefix}.prg")
    io_utils.write_prg(prefix, builder.prg)
    builder.serialize(f"{prefix}.pickle")
    # m = aseq.max_nesting_level_reached
    # logging.info(f"Max_nesting_reached\t{m}")

    # logging.info(f"Write GFA file to {prefix}.gfa")
    # io_utils.write_gfa(f"{prefix}.gfa", aseq.prg)



def run(cl_options):
    global options
    options = cl_options
    input_files = get_all_input_files(options.input)
    os.makedirs(options.output, exist_ok=True)
    with Pool(options.threads) as pool:
        pool.map(process_MSA, input_files, chunksize=1)
