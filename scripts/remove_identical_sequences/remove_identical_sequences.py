import argparse
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Remove identical sequences from fasta file.')
parser.add_argument('--fasta', type=str, help='Input fasta files', required=True)
args = parser.parse_args()

sequences_already_output = set()

with open(f"{args.fasta}.unique.fa", "w") as unique_sequences_fh, \
     open(f"{args.fasta}.repeated.fa", "w") as repeated_sequences_fh:
    for record in SeqIO.parse(args.fasta, "fasta"):
        if record.seq not in sequences_already_output:
            SeqIO.write(record, unique_sequences_fh, "fasta")
            sequences_already_output.add(record.seq)
        else:
            SeqIO.write(record, repeated_sequences_fh, "fasta")
