import logging
import os
from pathlib import Path
import shutil
import multiprocessing
import sys

from make_prg import io_utils
from make_prg.prg_builder import PrgBuilderCollection
from make_prg.utils import output_files_already_exist, setup_stderr_logging
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
    subparser_update_prg.add_argument(
        "--keep_temp",
        action="store_true",
        default=False,
        help="Keep temp files."
    )

    subparser_update_prg.set_defaults(func=run)

    return subparser_update_prg

def update(locus_name, variant_nodes_with_mutation, prg_builder_for_locus, temp_dir):
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
            except RuntimeError as exc:
                print("Failed finding leaf", file=sys.stderr)
                print(exc, file=sys.stderr)
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

    logging.info(f"Reading {options.pickle}...")
    prg_builder_collection = PrgBuilderCollection.deserialize(options.pickle)
    logging.info(f"Reading {options.denovo_paths}...")
    denovo_paths_db = DenovoPathsDB(options.denovo_paths)

    output_dir = Path(options.output_prefix).parent
    os.makedirs(output_dir, exist_ok=True)
    temp_path = Path(options.output_prefix+"_tmp")
    os.makedirs(temp_path, exist_ok=True)

    # update all PRGs with denovo sequences
    logging.info(f"Using {options.threads} threads to update PRGs...")
    multithreaded_input = []
    for locus_name, prg_builder_for_locus in prg_builder_collection.locus_name_to_prg_builder.items():  # we do for all PRGs as those that don't have denovo variants will be generated also
        variant_nodes_with_mutation = denovo_paths_db.locus_name_to_variant_nodes_with_mutation.get(locus_name, [])
        multithreaded_input.append((locus_name, variant_nodes_with_mutation, prg_builder_for_locus, temp_path))

    # avoids multiprocessing Pool deadlocks (see https://pythonspeed.com/articles/python-multiprocessing/)
    multiprocessing.set_start_method("spawn")
    with multiprocessing.Pool(options.threads) as pool:
        pool.starmap(update, multithreaded_input, chunksize=1)
    logging.info(f"All PRGs updated!")

    # concatenate output PRGs
    logging.info("Concatenating files from several threads into single final files...")
    # TODO: do this glob optimisation also in from msa
    prg_files = [f"{temp_path}/{locus_name}/{locus_name}.prg.fa" for locus_name in prg_builder_collection.locus_name_to_prg_builder.keys()]
    io_utils.concatenate_text_files(prg_files, options.output_prefix + ".prg.fa")

    stats_files = [f"{temp_path}/{locus_name}/{locus_name}.stats" for locus_name in
                 prg_builder_collection.locus_name_to_prg_builder.keys()]
    io_utils.concatenate_text_files(stats_files, options.output_prefix + ".stats")

    # remove temp files if needed
    if not options.keep_temp and temp_path.exists():
        logging.info("Removing temp files...")
        shutil.rmtree(temp_path)

    logging.info("All done!")
