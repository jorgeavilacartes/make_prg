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


sample_prg = " 5 AATAGGCCG 7  9 GATGCAGTTCAA 10 GATGCGGCGTA 9 AACGCCTTATCCGGCATACGA 11 ATTTATT 12 TTTTATT 11  8  13 G " \
             "14 A 13 ATGCGGCGTACGAATTTAT 15 T 16 C 15  7 CGGCCTGGCTCCCCGTAGGCCG 17 A 18 G 17 " \
             "ATAAGATGCGCCAGCATCGCATCCGGCTATAATGC 19 G 20 A 19  6  21 TTCATTGG 22 TTCAATG 23 G 24 A 23  21 " \
             "TTTATAATGCCTGATAAACGCACGGTCGATCCCCTCGCCCCTTCGGGGAGAGGATTAGGGTGAGGGGGTACAAGCCAGCCAGAGACCAGGCAA 25 " \
             "TGACATG 26 CGACATG 25  5 CACATAACC 27 TCT 28 ACC 27 TGAAACT 29  31 CTTT 32 CGTC 32 CTTC 32 CATC 31 " \
             "CCCAGAGCCTCTT 33 CAGC 34 TAGC 34 CAGT 33 CATCTATT 35 CA 36 TG 35  30  37 ACATCTCTTCA 38 ACATTTCTTCA 37  " \
             "29 GGAGCAAACAATTTCAT 39 GCCAACTC 40 TCCAACTC 40 TCCAACTT 40 ACCAACTC 39 ATAACCCCAGCATATAAATCCAG 41 T 42 " \
             "A 41 TGGTAACTTTT 43 A 44 C 43 TTTAACCT 45 G 46 A 45 AAACCAGTTT 47 TATCCAC 49 T 50 C 49  48 AATCCACC 47 " \
             "ATTTATAAAATTATGTGAAGCATTTCATAGAAGAAAAATCACTGGC 51 C 52 T 51 TAAACATTAT 53 C 54 T 53 CCCTTTTTGC 55 CTGG " \
             "56 CTGA 56 ATGA 56 CTAG 56 CTTA 56 CTGT 55 TTTTTGACCATTTCCG 57 C 58 T 57 " \
             "GATTTGTTACACATTGAAATATCACTTTTGCTGTGCGTAATATGGCTATTCGTTAGC 59 C 60 A 59 AAAAAATAAGAAAAGAT 61 T 62 A 61 "
