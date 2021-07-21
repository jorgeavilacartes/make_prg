########################################################################################################################
# configs
pandora_out = "pandora_out"
output_dir = "fastas"
########################################################################################################################

import gzip
import os
from glob import glob

os.makedirs(output_dir)
files = glob(f"{pandora_out}/**/pandora.consensus.fq.gz")
for fasta_file in files:
    sample = fasta_file.split("/")[-2]
    with gzip.open(fasta_file) as fasta_fh:
        while True:
            fasta_header, sequence, _, __ = fasta_fh.readline(), fasta_fh.readline(), fasta_fh.readline(), fasta_fh.readline()

            finished_reading = len(fasta_header) == 0
            if finished_reading:
                break

            fasta_header = fasta_header.decode("utf-8")
            sequence = sequence.decode("utf-8")
            locus_file_name = fasta_header.split()[0][1:]
            with open(f"{output_dir}/{locus_file_name}.fa", "a") as locus_fh:
                locus_fh.write(f">{sample}\n")
                locus_fh.write(sequence)