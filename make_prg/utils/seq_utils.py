from typing import Generator, Tuple, Iterable, List
import itertools
from Bio import pairwise2
from loguru import logger
from make_prg.from_msa import MSA
import copy
import numpy as np
from Bio.Seq import Seq


NONMATCH = "*"
GAP = "-"
Sequence = str
Sequences = List[str]


def is_non_match(letter: str):
    return letter == NONMATCH


def is_gap(letter: str):
    return letter == GAP


def remove_duplicates(seqs: Sequences) -> Generator:
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


def count(iterable) -> int:
    return sum(1 for _ in iterable)


def get_alignment_seqs(alignment: MSA) -> Generator:
    for record in alignment:
        yield str(record.seq)


def get_number_of_unique_ungapped_sequences(sub_alignment: MSA) -> int:
    alignment_seqs = list(get_alignment_seqs(sub_alignment))
    ungapped_sequences = list(map(ungap, alignment_seqs))
    deduplicated_ungapped_sequences = list(remove_duplicates(ungapped_sequences))
    number_of_unique_nongapped_sequences = count(deduplicated_ungapped_sequences)
    return number_of_unique_nongapped_sequences


def get_number_of_unique_gapped_sequences(sub_alignment: MSA) -> int:
    alignment_seqs = get_alignment_seqs(sub_alignment)
    deduplicated_gapped_sequences = remove_duplicates(alignment_seqs)
    number_of_unique_gapped_sequences = count(deduplicated_gapped_sequences)
    return number_of_unique_gapped_sequences


def get_expanded_sequences(alignment: MSA) -> Sequences:
    """
    Replace - with nothing, remove seqs containing N or other non-allowed letters
    and duplicate sequences containing RYKMSW, replacing with AGCT alternatives

    The sequences are deliberately returned in the order they are received
    """
    gapless_seqs = map(ungap, get_alignment_seqs(alignment))

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


def remove_gaps_from_MSA(alignment: MSA) -> MSA:
    """
    Return a gapless alignment. This code is long and a bit convoluted because it is optimised (it was too slow if
    done in the most intuitive way).
    """
    alignment_as_array = np.array([list(rec) for rec in alignment], str, order="F")
    gapless_sequences = [[] for _ in range(len(alignment))]
    for column_index in range(alignment.get_alignment_length()):
        column_bases = alignment_as_array[:, column_index]
        column_bases_deduplicated = list(set(column_bases))
        just_gaps = column_bases_deduplicated == [GAP]
        if not just_gaps:
            for gapless_sequence, base in zip(gapless_sequences, column_bases):
                gapless_sequence.append(base)

    gapless_records = []
    for gapless_sequence, previous_record in zip(gapless_sequences, alignment):
        new_record = copy.deepcopy(previous_record)
        new_record.seq = Seq("".join(gapless_sequence))
        gapless_records.append(new_record)

    gapless_alignment = MSA(gapless_records)
    return gapless_alignment


def get_consensus_from_MSA(alignment: MSA) -> str:
    """Produces a 'consensus string' from an MSA: at each position of the
    MSA, the string has a base if all aligned sequences agree, and a "*" if not.
    IUPAC ambiguous bases result in non-consensus and are later expanded in the prg.
    N results in consensus at that position unless they are all N."""
    consensus_string_as_list = []
    for i in range(alignment.get_alignment_length()):
        column = set([record.seq[i] for record in alignment])
        column = column.difference({"N"})
        if len(ambiguous_bases.intersection(column)) > 0 or len(column) != 1:
            consensus_string_as_list.append(NONMATCH)
        else:
            consensus_string_as_list.append(column.pop())
    consensus_string = "".join(consensus_string_as_list)
    return consensus_string
