from typing import List
from collections import namedtuple, defaultdict
import re

MLPathNode = namedtuple("MLPathNode", ["start_index", "end_index", "sequence"])
DenovoVariant = namedtuple("DenovoVariant", ["start_index", "ref", "alt"])


class DenovoLocusInfo:
    def __init__(self, sample: str, locus: str, ml_path: List[MLPathNode], variants: List[DenovoVariant]):
        self.sample = sample
        self.locus = locus
        self.ml_path = ml_path
        self.variants = variants


class DenovoPathsDB:
    def __read_ml_path(self, filehandler):
        line = filehandler.readline().strip()
        nb_of_nodes_in_ml_path = int(line.split()[0])
        ml_path = []
        for _ in range(nb_of_nodes_in_ml_path):
            line = filehandler.readline().strip()
            matches = self.ml_path_regex.search(line)
            ml_path.append(MLPathNode(
                start_index=matches.group(1),
                end_index=matches.group(2),
                sequence=matches.group(3)))
        return ml_path

    def __read_variants(self, filehandler):
        line = filehandler.readline().strip()
        nb_of_variants = int(line.split()[0])
        variants = []
        for _ in range(nb_of_variants):
            line = filehandler.readline().strip()
            line_split = line.split()
            variants.append(DenovoVariant(
                start_index=int(line_split[0]),
                ref=line_split[1],
                alt=line_split[2]
            ))
        return variants

    def __init__(self, filename):
        # Example:
        # (0 [0, 110) ATGCAGATACGTGAACAGGGCCGCAAAATTCAGTGCATCCGCACCGTGTACGACAAGGCCATTGGCCGGGGTCGGCAGACGGTCATTGCCACACTGGCCCGCTATACGAC)
        self.ml_path_regex = re.compile("\(\d+ \[(\d+), (\d+)\) ([ACGT]+)\)")

        self._locus_name_to_denovo_loci = defaultdict(list)
        with open(filename) as filehandler:
            # read nb_of_samples
            line = filehandler.readline().strip()
            nb_of_samples = int(line.split()[0])

            for sample_index in range(nb_of_samples):
                # read each sample
                line = filehandler.readline().strip()
                sample = line.split()[1]
                line = filehandler.readline().strip()
                nb_of_loci_in_sample = int(line.split()[0])
                for locus_index in range(nb_of_loci_in_sample):
                    # read each locus
                    locus = filehandler.readline().strip()
                    ml_path = self.__read_ml_path(filehandler)
                    variants = self.__read_variants(filehandler)
                    denovo_locus = DenovoLocusInfo(sample, locus, ml_path, variants)
                    self._locus_name_to_denovo_loci[locus].append(denovo_locus)

