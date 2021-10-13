from typing import Generator, Sequence, Tuple
import itertools
from Bio import pairwise2

from loguru import logger

from make_prg.from_msa import MSA

NONMATCH = "*"
GAP = "-"


def is_non_match(letter: str):
    return letter == NONMATCH


def is_gap(letter: str):
    return letter == GAP


def remove_duplicates(seqs: Sequence) -> Generator:
    seen = set()
    for x in seqs:
        if x in seen:
            continue
        seen.add(x)
        yield x


iupac = {
    "R": "GA",
    "Y": "TC",
    "K": "GT",
    "M": "AC",
    "S": "GC",
    "W": "AT",
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
}
allowed_bases = set(iupac.keys())
standard_bases = {"A", "C", "G", "T"}
ambiguous_bases = allowed_bases.difference(standard_bases)


def has_empty_sequence(alignment: MSA, interval: Tuple[int, int]) -> bool:
    sub_alignment = alignment[:, interval[0] : interval[1] + 1]
    for record in sub_alignment:
        if all(map(is_gap, record.seq)):
            return True
    return False


def ungap(seq: str) -> str:
    return seq.replace("-", "")


def get_interval_seqs(interval_alignment: MSA):
    """
    Replace - with nothing, remove seqs containing N or other non-allowed letters
    and duplicate sequences containing RYKMSW, replacing with AGCT alternatives

    The sequences are deliberately returned in the order they are received
    """
    gapless_seqs = [ungap(str(record.seq)) for record in interval_alignment]

    callback_seqs, expanded_seqs = [], []
    expanded_set = set()
    for seq in remove_duplicates(gapless_seqs):
        if len(expanded_set) == 0:
            callback_seqs.append(seq)
        if not set(seq).issubset(allowed_bases):
            continue
        alternatives = [iupac[base] for base in seq]
        for tuple_product in itertools.product(*alternatives):
            expanded_str = "".join(tuple_product)
            if expanded_str not in expanded_set:
                expanded_set.add(expanded_str)
                expanded_seqs.append(expanded_str)

    if len(expanded_set) == 0:
        logger.warning(
            "Every sequence must have contained an N in this slice - redo sequence curation"
        )
        logger.warning(f'Sequences were: {" ".join(callback_seqs)}')
        logger.warning("Using these sequences anyway, and should be ignored downstream")
        return callback_seqs
    return expanded_seqs


def align(ref: str, alt: str, match_score=2, mismatch_score=-1, gap_open_score=-4, gap_extend_score=-2) \
        -> Tuple[str, str]:
    alignment = pairwise2.align.globalms(
        ref, alt, match_score, mismatch_score, gap_open_score, gap_extend_score
    )
    empty_alignment = len(alignment) == 0
    if empty_alignment:
        #  this usually happens if ref or alt are empty, let's check
        ref_is_empty = len(ref) == 0
        alt_is_empty = len(alt) == 0

        # only one should be empty
        both_are_empty = ref_is_empty and alt_is_empty
        assert not both_are_empty
        both_are_not_empty = (not ref_is_empty) and (not alt_is_empty)
        assert not both_are_not_empty

        if ref_is_empty:
            return GAP * len(alt), alt
        elif alt_is_empty:
            return ref, GAP * len(ref)
        else:
            assert True, "Unreachable code"  # just to be sure
    else:
        alignment = alignment[
            0
        ]  # get the first alignment - we might have several equally good ones
        return alignment.seqA, alignment.seqB
