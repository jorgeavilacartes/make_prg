from pathlib import Path
from collections import defaultdict
import pysam

pandora_discovery_dir = Path(snakemake.input.pandora_discovery_dir).absolute()
sample_dirs = [sample for sample in pandora_discovery_dir.iterdir()
               if sample.is_dir()]

locus_to_sample_to_seq = defaultdict(lambda: defaultdict(str))
for sample_dir in sample_dirs:
    sample = sample_dir.name
    with pysam.FastxFile(str(sample_dir / "denovo_sequences.fa")) as denovo_seq_fh:
        for fasta_record in denovo_seq_fh:
            locus = fasta_record.name
            sequence = fasta_record.sequence
            locus_to_sample_to_seq[locus][sample] = sequence

fasta_output_dir = Path(f'{snakemake.config["output_dir"]}/fastas/')
fasta_output_dir.mkdir(parents=True)
for locus, sample_to_seq in locus_to_sample_to_seq.items():
    with open(f"{fasta_output_dir}/{locus}.fa", "w") as locus_fh:
        for sample, seq in sample_to_seq.items():
            print(f">{locus}_{sample}", file=locus_fh)
            print(seq, file=locus_fh)
