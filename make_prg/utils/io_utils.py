import re
import gzip
from pathlib import Path
import fileinput
import multiprocessing
from Bio import AlignIO
import os
from loguru import logger
from make_prg.from_msa import MSA
from typing import Dict, Optional, List, Tuple
from collections import defaultdict
from make_prg import prg_builder
from make_prg.utils.prg_encoder import PrgEncoder, PRG_Ints
from pathlib import Path
from zipfile import ZipFile


def load_alignment_file(msa_file: str, alignment_format: str) -> MSA:
    msa_file = str(msa_file)
    if msa_file.endswith(".gz"):
        handle = gzip.open(msa_file, "rt")
        alignment = AlignIO.read(handle, alignment_format)
        handle.close()
    else:
        alignment = AlignIO.read(msa_file, alignment_format)
    for record in alignment:
        record.seq = record.seq.upper()
    return alignment


# ************/
# GFA code */
# ***********/
class GFA_Output:
    """
    A simple class for converting a PRG string into a GFA file
    TODO: Update to GFA2 format
    """

    def __init__(self, gfa_string="", gfa_id=0, gfa_site=5):
        self.gfa_string = gfa_string
        self.gfa_id = gfa_id
        self.gfa_site = gfa_site
        self.delim_char = " "  # This mirrors the AlignedSeq class.

    def split_on_site(self, prg_string, site_num):
        site_coords = [
            (a.start(), a.end())
            for a in list(
                re.finditer(
                    "%s%d%s" % (self.delim_char, site_num, self.delim_char), prg_string
                )
            )
        ]
        last_pos = None
        split_strings = []
        for (start, end) in site_coords:
            split_strings.append(prg_string[last_pos:start])
            last_pos = end
        split_strings.append(prg_string[last_pos:])
        delim = "%s%d%s" % (self.delim_char, site_num, self.delim_char)
        check_string = delim.join(split_strings)
        assert check_string == prg_string, (
                "Something has gone wrong with the string split for site %d\nsplit_"
                "strings: %s" % (site_num, split_strings)
        )
        return split_strings

    def build_gfa_string(self, prg_string, pre_var_id=None):
        """Takes prg_string and builds a gfa_string with fragments
        from the prg_string."""
        end_ids = []
        # iterate through sites present, updating gfa_string with each in turn
        while str(self.gfa_site) in prg_string:
            logger.trace("gfa_site: {}", self.gfa_site)
            prgs = self.split_on_site(prg_string, self.gfa_site)
            logger.trace("prgs: {}", prgs)
            assert len(prgs) == 3, "Invalid prg sequence %s for site %d and id %d" % (
                prg_string,
                self.gfa_site,
                self.gfa_id,
            )

            # add pre-var site string and links from previous seq fragments
            if prgs[0] != "":
                self.gfa_string += "S\t%d\t%s\tRC:i:0\n" % (self.gfa_id, prgs[0])
            else:
                # adds an empty node for empty pre var site seqs
                self.gfa_string += "S\t%d\t%s\tRC:i:0\n" % (self.gfa_id, "*")
            pre_var_id = self.gfa_id
            self.gfa_id += 1
            for id in end_ids:
                self.gfa_string += "L\t%d\t+\t%d\t+\t0M\n" % (id, pre_var_id)
                end_ids = []

            # recursively add segments for each of the variant haplotypes at
            # this site, saving the end id for each haplotype
            vars = self.split_on_site(prgs[1], self.gfa_site + 1)
            assert len(vars) > 1, "Invalid prg sequence %s for site %d and id %d" % (
                prg_string,
                self.gfa_site + 1,
                self.gfa_id,
            )
            logger.trace("vars: {}", vars)
            self.gfa_site += 2
            logger.trace("gfa_site: {}", self.gfa_site)
            for var_string in vars:
                if pre_var_id != None:
                    self.gfa_string += "L\t%d\t+\t%d\t+\t0M\n" % (
                        pre_var_id,
                        self.gfa_id,
                    )
                var_end_ids = self.build_gfa_string(
                    prg_string=var_string, pre_var_id=pre_var_id
                )
                end_ids.extend(var_end_ids)

            prg_string = prgs[2]
            pre_var_id = None

        # finally add the final bit of sequence after variant site
        if prg_string != "":
            self.gfa_string += "S\t%d\t%s\tRC:i:0\n" % (self.gfa_id, prg_string)
        else:
            self.gfa_string += "S\t%d\t%s\tRC:i:0\n" % (self.gfa_id, "*")
        for id in end_ids:
            self.gfa_string += "L\t%d\t+\t%d\t+\t0M\n" % (id, self.gfa_id)
        end_ids = []
        return_id = [self.gfa_id]
        self.gfa_id += 1
        return return_id


def write_gfa(prefix, prg_string):
    """
    Writes a gfa file from prg string.
    """
    with open(f"{prefix}.gfa", "w") as f:
        # initialize gfa_string, id and site, then update string with the prg
        gfa_string = "H\tVN:Z:1.0\tbn:Z:--linear --singlearr\n"
        gfa_obj = GFA_Output(gfa_string)
        gfa_obj.build_gfa_string(prg_string=prg_string)
        f.write(gfa_obj.gfa_string)


# ******************/
# Write PRG code */
# *****************/


def write_prg(output_prefix: str, prg_string: str):
    """
    Writes the prg to outfile.
    Writes it as a human readable string, and also as an integer vector
    """
    sample = Path(output_prefix).with_suffix("").name
    prg_filename = Path(output_prefix + ".prg.fa")
    with prg_filename.open("w") as prg:
        print(f">{sample}\n{prg_string}", file=prg)

    prg_ints_fpath = Path(output_prefix + ".bin")
    prg_encoder = PrgEncoder()
    prg_ints: PRG_Ints = prg_encoder.encode(prg_string)
    with prg_ints_fpath.open("wb") as ostream:
        prg_encoder.write(prg_ints, ostream)


def concatenate_text_files(input_filepaths, output_filepath):
    output_filepath_parent_dir = Path(output_filepath).parent
    os.makedirs(output_filepath_parent_dir, exist_ok=True)

    with open(output_filepath, "w") as fout:
        empty_input = len(input_filepaths) == 0
        if empty_input:
            return
        else:
            with fileinput.input(input_filepaths) as fin:
                for line in fin:
                    fout.write(line)


# From https://gist.github.com/jacobtomlinson/9031697
def remove_empty_folders(path: str, remove_root: bool = True):
    if not os.path.isdir(path):
        return

    # remove empty subfolders
    files = os.listdir(path)
    if len(files):
        for f in files:
            fullpath = os.path.join(path, f)
            if os.path.isdir(fullpath):
                remove_empty_folders(fullpath)

    # if folder empty, delete it
    files = os.listdir(path)
    if len(files) == 0 and remove_root:
        os.rmdir(path)


def output_files_already_exist(output_prefix: str):
    return (
        Path(output_prefix + ".prg.fa").exists()
        or Path(output_prefix + ".update_DS.zip").exists()
        or Path(output_prefix + ".prg.bin.zip").exists()
    )


def get_temp_dir_for_multiprocess(root_temp_dir: Path):
    current_process = multiprocessing.current_process()
    temp_dir = root_temp_dir / "mp_temp" / current_process.name
    os.makedirs(temp_dir, exist_ok=True)
    return temp_dir


# get all files that were generated
class SetOutputFiles:
    def __init__(self, PRG: Optional[Path] = None, binary_PRG: Optional[Path] = None,
                 gfa: Optional[Path] = None, pickle: Optional[Path] = None,
                 stats: Optional[Path] = None):
        self.PRG: Optional[Path] = PRG
        self.binary_PRG: Optional[Path] = binary_PRG
        self.gfa: Optional[Path] = gfa
        self.pickle: Optional[Path] = pickle
        self.stats: Optional[Path] = stats

    @staticmethod
    def clear_temp_extensions(filename: str) -> str:
        while True:
            file_should_be_cleared = filename.endswith(".fa") or filename.endswith(".prg") or \
                                     filename.endswith(".pickle") or filename.endswith(".stats") or \
                                     filename.endswith(".bin") or filename.endswith(".gfa")
            if file_should_be_cleared:
                filename = Path(filename).with_suffix("").name
            else:
                return filename

    @staticmethod
    def _delete_file(filepath: Optional[Path]):
        if filepath is not None:
            filepath.unlink(missing_ok=True)

    def delete_files(self):
        self._delete_file(self.PRG)
        self._delete_file(self.binary_PRG)
        self._delete_file(self.gfa)
        self._delete_file(self.pickle)
        self._delete_file(self.stats)


def get_locus_to_set_of_output_files(threads: int, temp_root: Path) -> Dict[str, SetOutputFiles]:
    locus_to_set_of_output_files = defaultdict(SetOutputFiles)
    for process_num in range(1, threads + 1):
        workdir = temp_root / f"ForkPoolWorker-{process_num}"
        if workdir.exists():
            for file in workdir.iterdir():
                if file.is_file():
                    locus_name = SetOutputFiles.clear_temp_extensions(file.name)
                    if file.name.endswith(".prg.fa"):
                        locus_to_set_of_output_files[locus_name].PRG = file
                    elif file.name.endswith(".bin"):
                        locus_to_set_of_output_files[locus_name].binary_PRG = file
                    elif file.name.endswith(".gfa"):
                        locus_to_set_of_output_files[locus_name].gfa = file
                    elif file.name.endswith(".pickle"):
                        locus_to_set_of_output_files[locus_name].pickle = file
                    elif file.name.endswith(".stats"):
                        locus_to_set_of_output_files[locus_name].stats = file
    return locus_to_set_of_output_files


def get_stats_on_variants(stats_files: List[Path]) -> Tuple[int, int]:
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


def zip_set_of_files(zip_filepath: Path, filename_to_filepath: Dict[str, Path]):
    is_a_zip_file = zip_filepath.suffix == ".zip"
    assert is_a_zip_file, "zip_set_of_files() was not given a .zip filepath"
    with ZipFile(zip_filepath, "w") as zip_file:
        for filename, filepath in filename_to_filepath.items():
            zip_file.write(filepath, filename)


def create_final_files(threads: int, output_prefix: str, output_stats: bool = False):
    logger.info("Concatenating files from several threads into single final files...")

    logger.info("Creating FASTA file of PRGs...")
    locus_to_set_of_output_files = get_locus_to_set_of_output_files(threads, Path(output_prefix) / "mp_temp")

    prg_files = [output_files.PRG for output_files in locus_to_set_of_output_files.values()]
    concatenate_text_files(prg_files, output_prefix + ".prg.fa")

    # zip all PRG Builders
    logger.info("Creating zip of update data structures...")
    prg_builder_zip_db = prg_builder.PrgBuilderZipDatabase(Path(f"{output_prefix}.update_DS.zip"))
    locus_to_prg_builder_pickle_path = {locus: output_files.pickle
                                        for locus, output_files in locus_to_set_of_output_files.items()}
    prg_builder_zip_db.save(locus_to_prg_builder_pickle_path)

    # zip all encoded PRGs
    logger.info("Creating zip of encoded PRGs...")
    filename_to_encoded_PRG_paths = {output_files.binary_PRG.name: output_files.binary_PRG
                                     for output_files in locus_to_set_of_output_files.values()}
    zip_set_of_files(Path(f"{output_prefix}.prg.bin.zip"), filename_to_encoded_PRG_paths)

    # zip all GFAs
    logger.info("Creating zip of GFAs...")
    filename_to_gfa_paths = {output_files.gfa.name: output_files.gfa
                             for output_files in locus_to_set_of_output_files.values()}
    zip_set_of_files(Path(f"{output_prefix}.prg.gfa.zip"), filename_to_gfa_paths)

    # sum up stats files and output stats
    if output_stats:
        logger.info("Computing stats on updates...")
        stats_files = [output_files.stats for output_files in locus_to_set_of_output_files.values()]
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

    # cleanup
    for output_files in locus_to_set_of_output_files.values():
        output_files.delete_files()
    remove_empty_folders(output_prefix)
