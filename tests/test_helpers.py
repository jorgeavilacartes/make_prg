from typing import List
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from make_prg.from_msa import MSA


def make_alignment(seqs: List[str], ids: List[str] = None) -> MSA:
    seq_lengths = set(map(len, seqs))
    assert (
        len(seq_lengths) == 1
    ), "Sequences are not the same length, does not represent an alignment"
    if ids is None:
        seqrecords = [SeqRecord(Seq(seq), id=f"s{i}", description="") for i, seq in enumerate(seqs)]
    else:
        seqrecords = [SeqRecord(Seq(seq), id=ID, description="") for seq, ID in zip(seqs, ids)]
    return MSA(seqrecords)


def equal_msas(msa_1: MSA, msa_2: MSA) -> bool:
    msa_1_as_fasta = format(msa_1, "fasta")
    msa_2_as_fasta = format(msa_2, "fasta")
    return msa_1_as_fasta == msa_2_as_fasta
