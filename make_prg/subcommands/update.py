import logging
import os
from pathlib import Path
import shutil
import multiprocessing
import sys

from make_prg import io_utils
from make_prg.prg_builder import PrgBuilderCollection, PrgBuilder, LeafNotFoundException
from make_prg.utils import output_files_already_exist, setup_stderr_logging
from make_prg.denovo_paths_reader import DenovoPathsDB

def register_parser(subparsers):
    subparser_update_prg = subparsers.add_parser(
        "update",
        usage="make_prg update_prg",
        help="Update PRGs given new sequences output by pandora.",
    )
    subparser_update_prg.add_argument(
        "-i", "--input_prefix",
        action="store",
        type=str,
        required=True,
        help=(
            "Input filepath prefix to update data structures. Points to <prefix>.update_DS.bak, <prefix>.update_DS.dat, and "
            "<prefix>.update_DS.dir."
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
    subparser_update_prg.add_argument(
        "--keep_temp",
        action="store_true",
        default=False,
        help="Keep temp files."
    )

    subparser_update_prg.set_defaults(func=run)

    return subparser_update_prg

def get_stats_on_variants(stats_files):
    nb_of_variants_successfully_applied = 0
    nb_of_variants_that_failed_to_be_applied = 0
    for stat_file in stats_files:
        with open(stat_file) as stat_file_fh:
            line_split = stat_file_fh.readline().strip().split()
            nb_of_variants_successfully_applied_for_this_locus = int(line_split[1])
            nb_of_variants_successfully_applied += nb_of_variants_successfully_applied_for_this_locus
            nb_of_variants_that_failed_to_be_applied_for_this_locus = int(line_split[2])
            nb_of_variants_that_failed_to_be_applied += nb_of_variants_that_failed_to_be_applied_for_this_locus
    return nb_of_variants_successfully_applied, nb_of_variants_that_failed_to_be_applied


def update(locus_name, prg_builder_pickle_filepath, variant_nodes_with_mutation, temp_dir):
    prg_builder_for_locus = PrgBuilder.deserialize(prg_builder_pickle_filepath)
    nb_of_variants_sucessfully_updated = 0
    nb_of_variants_with_failed_update = 0

    we_have_variants = len(variant_nodes_with_mutation) > 0
    if we_have_variants:
        print(f"Updating {locus_name} ...")

        leaves_to_update = set()
        for variant_node_with_mutation in variant_nodes_with_mutation:
            try:
                prg_builder_tree_node = prg_builder_for_locus.get_node_given_interval(variant_node_with_mutation.key)
                prg_builder_tree_node.add_seq_to_batch_update(variant_node_with_mutation.mutated_node_sequence)
                leaves_to_update.add(prg_builder_tree_node)
                nb_of_variants_sucessfully_updated += 1
            except LeafNotFoundException as exc:
                print(f"Failed finding leaf: {exc}", file=sys.stderr)
                nb_of_variants_with_failed_update += 1

        # update the changed leaves
        for leaf in leaves_to_update:
            leaf.batch_update(temp_dir)
        print(f"Updated {locus_name}: {len(variant_nodes_with_mutation)} denovo sequences added!")
    else:
        print(f"{locus_name} has no new variants, no update needed")

    # regenerate PRG
    locus_prefix = temp_dir / locus_name / locus_name
    locus_prefix_parent = locus_prefix.parent
    os.makedirs(locus_prefix_parent, exist_ok=True)
    prg = prg_builder_for_locus.build_prg()
    print(f"Write PRG file to {locus_prefix}.prg.fa")
    io_utils.write_prg(str(locus_prefix), prg)

    with open(f"{locus_prefix}.stats", "w") as stats_filehandler:
        print(f"{locus_name} {nb_of_variants_sucessfully_updated} {nb_of_variants_with_failed_update}", file=stats_filehandler)

    # Note: we intentionally do not regenerate updateable data structure here because we don't want to update
    # PRGs on top of already updated PRGs
    # TODO: change this?


def run(options):
    if output_files_already_exist(options.output_prefix):
        raise RuntimeError("One or more output files already exists, aborting run...")

    setup_stderr_logging()

    logging.info(f"Reading update data structures...")
    prg_builder_collection = PrgBuilderCollection.deserialize(f"{options.input_prefix}.update_DS")
    prg_builder_collection.to_absolute_paths()
    logging.info(f"Reading {options.denovo_paths}...")
    denovo_paths_db = DenovoPathsDB(options.denovo_paths)

    output_dir = Path(options.output_prefix).parent
    os.makedirs(output_dir, exist_ok=True)
    temp_path = Path(options.output_prefix+"_tmp")
    os.makedirs(temp_path, exist_ok=True)

    # update all PRGs with denovo sequences
    logging.info(f"Using {options.threads} threads to update PRGs...")
    multithreaded_input = []
    for locus_name, prg_builder_pickle_filepath in prg_builder_collection.locus_name_to_pickle_files.items():  # we do for all PRGs as those that don't have denovo variants will be generated also
        variant_nodes_with_mutation = denovo_paths_db.locus_name_to_variant_nodes_with_mutation.get(locus_name, [])
        multithreaded_input.append((locus_name, prg_builder_pickle_filepath, variant_nodes_with_mutation, temp_path))

    # avoids multiprocessing Pool deadlocks (see https://pythonspeed.com/articles/python-multiprocessing/)
    multiprocessing.set_start_method("spawn")
    with multiprocessing.Pool(options.threads) as pool:
        pool.starmap(update, multithreaded_input, chunksize=1)
    logging.info(f"All PRGs updated!")

    # concatenate output PRGs
    logging.info("Concatenating files from several threads into single final files...")
    prg_files = [f"{temp_path}/{locus_name}/{locus_name}.prg.fa" for locus_name in prg_builder_collection.locus_name_to_pickle_files.keys()]
    io_utils.concatenate_text_files(prg_files, options.output_prefix + ".prg.fa")

    # sum up stats files and output stats
    stats_files = [f"{temp_path}/{locus_name}/{locus_name}.stats" for locus_name in
                 prg_builder_collection.locus_name_to_pickle_files.keys()]
    nb_of_variants_successfully_applied, nb_of_variants_that_failed_to_be_applied = \
        get_stats_on_variants(stats_files)
    print(f"Number of variants successfully applied: {nb_of_variants_successfully_applied}")
    print(f"Number of variants that failed to be applied: {nb_of_variants_that_failed_to_be_applied}")

    # remove temp files if needed
    if not options.keep_temp and temp_path.exists():
        logging.info("Removing temp files...")
        shutil.rmtree(temp_path)

    logging.info("All done!")
