from typing import List, Tuple
import multiprocessing
import os
import shutil
from pathlib import Path
from loguru import logger
from make_prg.utils import io_utils
from make_prg.update.denovo_variants import DenovoVariantsDB, UpdateData
from make_prg.prg_builder import PrgBuilderCollection, PrgBuilder, LeafNotFoundException
from make_prg.utils.utils import output_files_already_exist
from make_prg.utils.msa_aligner import MAFFT, MSAAligner


def register_parser(subparsers):
    subparser_update_prg = subparsers.add_parser(
        "update",
        usage="make_prg update",
        help="Update PRGs given new sequences.",
    )
    subparser_update_prg.add_argument(
        "-u",
        "--update_DS",
        action="store",
        type=str,
        required=True,
        help=(
            "Filepath to the update data structures. Should point to a file *.update_DS."
        ),
    )
    subparser_update_prg.add_argument(
        "-d",
        "--denovo_paths",
        action="store",
        type=str,
        required=True,
        help=(
            "Filepath containing denovo sequences. Should point to a denovo_paths.txt file."
        ),
    )
    subparser_update_prg.add_argument(
        "-o",
        "--output_prefix",
        action="store",
        type=str,
        required=True,
        help="Output prefix: prefix for the output files",
    )
    subparser_update_prg.add_argument(
        "-t",
        "--threads",
        action="store",
        type=int,
        default=1,
        help="Number of threads",
    )
    subparser_update_prg.add_argument(
        "--mafft",
        help="Path to MAFFT executable. By default, it is assumed to be on $PATH",
        default="mafft",
    )
    subparser_update_prg.add_argument(
        "--keep_temp", action="store_true", default=False, help="Keep temp files."
    )

    subparser_update_prg.set_defaults(func=run)

    return subparser_update_prg


def get_stats_on_variants(stats_files: List[str]) -> Tuple[int, int]:
    nb_of_variants_successfully_applied = 0
    nb_of_variants_that_failed_to_be_applied = 0
    for stat_file in stats_files:
        with open(stat_file) as stat_file_fh:
            line_split = stat_file_fh.readline().strip().split()
            nb_of_variants_successfully_applied_for_this_locus = int(line_split[1])
            nb_of_variants_successfully_applied += (
                nb_of_variants_successfully_applied_for_this_locus
            )
            nb_of_variants_that_failed_to_be_applied_for_this_locus = int(line_split[2])
            nb_of_variants_that_failed_to_be_applied += (
                nb_of_variants_that_failed_to_be_applied_for_this_locus
            )
    return nb_of_variants_successfully_applied, nb_of_variants_that_failed_to_be_applied


def update(
    locus_name: str,
    prg_builder_pickle_filepath: Path,
    update_data_list: List[UpdateData],
    msa_aligner: MSAAligner,
    temp_dir: Path
):
    prg_builder_for_locus = PrgBuilder.deserialize(prg_builder_pickle_filepath)
    prg_builder_for_locus.aligner = msa_aligner

    nb_of_variants_sucessfully_updated = 0
    nb_of_variants_with_failed_update = 0

    we_have_variants = len(update_data_list) > 0
    if we_have_variants:
        logger.debug(f"Updating {locus_name} ...")

        leaves_to_update = set()
        for update_data in update_data_list:
            try:
                prg_builder_tree_node = prg_builder_for_locus.get_node_given_interval(
                    update_data.ml_path_node_key
                )
                prg_builder_tree_node.add_seq_to_batch_update(
                    update_data.new_node_sequence
                )
                leaves_to_update.add(prg_builder_tree_node)
                nb_of_variants_sucessfully_updated += 1
            except LeafNotFoundException as exc:
                logger.warning(f"Failed finding leaf: {exc}")
                nb_of_variants_with_failed_update += 1

        # update the changed leaves
        for leaf in leaves_to_update:
            leaf.batch_update()
        logger.debug(
            f"Updated {locus_name}: {nb_of_variants_sucessfully_updated} denovo sequences added!"
        )
    else:
        logger.debug(f"{locus_name} has no new variants, no update needed")

    # regenerate PRG
    locus_prefix = temp_dir / locus_name / locus_name
    locus_prefix_parent = locus_prefix.parent
    os.makedirs(locus_prefix_parent, exist_ok=True)

    logger.info(f"Write PRG file to {locus_prefix}.prg.fa")
    prg = prg_builder_for_locus.build_prg()
    io_utils.write_prg(str(locus_prefix), prg)
    with open(f"{locus_prefix}.stats", "w") as stats_filehandler:
        print(
            f"{locus_name} {nb_of_variants_sucessfully_updated} {nb_of_variants_with_failed_update}",
            file=stats_filehandler,
        )

    # Note: we intentionally do not regenerate updateable data structure here because we don't want to update
    # PRGs on top of already updated PRGs
    # TODO: change this?


def run(options):
    if output_files_already_exist(options.output_prefix):
        raise RuntimeError("One or more output files already exists, aborting run...")

    # read input data
    logger.info("Checking Multiple Sequence Aligner...")
    temp_path = Path(options.output_prefix + "_tmp")
    os.makedirs(temp_path, exist_ok=True)
    mafft_aligner = MAFFT(executable=options.mafft, tmpdir=temp_path)
    logger.info("Reading update data structures...")
    prg_builder_collection = PrgBuilderCollection.deserialize(options.update_DS)
    logger.info(f"Reading {options.denovo_paths}...")
    denovo_variants_db = DenovoVariantsDB(options.denovo_paths)

    output_dir = Path(options.output_prefix).parent
    os.makedirs(output_dir, exist_ok=True)

    # update all PRGs with denovo sequences
    logger.info(f"Using {options.threads} threads to update PRGs...")
    multithreaded_input = []
    for (
        locus_name,
        prg_builder_pickle_filepath,
    ) in (
        prg_builder_collection.locus_name_to_pickle_filepaths.items()
    ):  # we do for all PRGs as those that don't have denovo variants will be generated also
        update_data = denovo_variants_db.locus_name_to_update_data.get(locus_name, [])
        multithreaded_input.append(
            (
                locus_name,
                Path(prg_builder_pickle_filepath),
                update_data,
                mafft_aligner,
                temp_path
            )
        )

    with multiprocessing.Pool(options.threads, maxtasksperchild=1) as pool:
        pool.starmap(update, multithreaded_input, chunksize=1)
    logger.success(f"All PRGs updated!")

    # concatenate output PRGs
    logger.info("Concatenating files from several threads into a single file...")
    prg_files = [
        f"{temp_path}/{locus_name}/{locus_name}.prg.fa"
        for locus_name in prg_builder_collection.locus_name_to_pickle_filepaths.keys()
    ]
    io_utils.concatenate_text_files(prg_files, options.output_prefix + ".prg.fa")

    # sum up stats files and output stats
    stats_files = [
        f"{temp_path}/{locus_name}/{locus_name}.stats"
        for locus_name in prg_builder_collection.locus_name_to_pickle_filepaths.keys()
    ]
    (
        nb_of_variants_successfully_applied,
        nb_of_variants_that_failed_to_be_applied,
    ) = get_stats_on_variants(stats_files)
    logger.success(
        f"Number of variants successfully applied: {nb_of_variants_successfully_applied}"
    )
    logger.warning(
        f"Number of variants that failed to be applied: {nb_of_variants_that_failed_to_be_applied}"
    )

    # remove temp files if needed
    if not options.keep_temp and temp_path.exists():
        logger.info("Removing temp files...")
        shutil.rmtree(temp_path)

    logger.success("All done!")
