from collections import Counter, defaultdict
from typing import List, Dict, Iterator, Union, Optional
from itertools import starmap, repeat, chain
from loguru import logger
import numpy as np
from sklearn.cluster import KMeans
from make_prg.from_msa import MSA
from make_prg.utils.misc import flatten_list
from make_prg.utils.seq_utils import ungap, Sequence, Sequences, SequenceExpander
from dataclasses import dataclass

IDs = List[str]
SeqToIDs = Dict[Sequence, IDs]
SeqToSeqs = Dict[Sequence, Sequences]
ClusteredIDs = List[IDs]
ClusteredSeqs = List[Sequences]
KmerIDs = Dict[Sequence, int]

DISTANCE_THRESHOLD: float = 0.2
LENGTH_THRESHOLD: int = 5
MAX_CLUSTERS: int = 10


def count_distinct_kmers(seqs: Sequences, kmer_size: int) -> KmerIDs:
    for seq in seqs:
        if len(seq) < kmer_size:
            raise ValueError(f"Input sequence {seq} has length < kmer size {kmer_size}")
    result = dict()
    kmer_id = 0
    for seq in seqs:
        for start_pos in range(len(seq) - kmer_size + 1):
            kmer = seq[start_pos : start_pos + kmer_size]
            if kmer not in result:
                result[kmer] = kmer_id
                kmer_id += 1
    return result


def count_kmer_occurrences(seqs: Sequences, kmers: KmerIDs):
    """
    Computes a count matrix of kmers in :param seqs
    :param kmers: all possible kmers assumed to have been counted in seqs
    :return The counts of all kmers in each input sequence
    """
    num_kmers = len(kmers)
    kmer_size = len(next(iter(kmers)))
    result = np.zeros(shape=(len(seqs), num_kmers))
    for j, seq in enumerate(seqs):
        counts = result[j]
        for i in range(len(seq) - kmer_size + 1):
            kmer = seq[i : i + kmer_size]
            counts[kmers[kmer]] += 1
        result[j] = counts
    return result


def get_majority_char_in_column(sequences: Sequences, col_idx: int) -> str:
    max_idx = len(sequences[0]) - 1
    if not (0 <= col_idx <= max_idx):
        raise ValueError(f"Column index {col_idx} not in range(0,{max_idx})")

    column = [seq[col_idx] for seq in sequences]
    return Counter(column).most_common(1)[0][0]


def get_majority_string(sequences: Sequences) -> str:
    if len(sequences) == 1:
        return sequences[0]
    seqlen = len(sequences[0])
    if not all([len(seq) == seqlen for seq in sequences]):
        raise ValueError("Not all sequences have the same length")

    arguments = zip(repeat(sequences), range(seqlen))
    return "".join(starmap(get_majority_char_in_column, arguments))


def hamming_distance(seq1: Sequence, seq2: Sequence) -> int:
    seqlen = len(seq1)
    result = 0
    for i in range(seqlen):
        if seq1[i] != seq2[i]:
            result += 1
    return result


def get_distances(sequences: Sequences, consensus_string: Sequence) -> Iterator[int]:
    arguments = zip(sequences, repeat(consensus_string))
    return starmap(hamming_distance, arguments)


def get_one_ref_like_threshold_distance(seqlen: int) -> int:
    if seqlen < LENGTH_THRESHOLD:
        return 1
    else:
        return int(DISTANCE_THRESHOLD * seqlen)


def sequences_are_one_reference_like(sequences: Sequences) -> bool:
    majority_string = get_majority_string(sequences)
    distance_iterator = get_distances(sequences, majority_string)
    min_thresh = get_one_ref_like_threshold_distance(len(sequences[0]))
    return all(map(lambda dist: dist <= min_thresh, distance_iterator))


def cluster_further(clusters: ClusteredSeqs) -> bool:
    not_ref_like = lambda sequences: not sequences_are_one_reference_like(sequences)
    return any(map(not_ref_like, clusters))


def extract_clusters(
    seqdict: Union[SeqToIDs, SeqToSeqs], cluster_assignment: List[int]
) -> Union[ClusteredIDs, ClusteredSeqs]:
    value_pool = list(seqdict.values())
    num_elems = len(cluster_assignment)
    if num_elems != len(value_pool):
        raise ValueError(
            f"Mismatch between number of sequences/ID lists and number of cluster assignments"
        )
    num_clusters = max(cluster_assignment) + 1
    if set(cluster_assignment) != set(range(num_clusters)):
        raise ValueError(
            "Inconsistent cluster numbering (likely reason: more input sequences that clustered data points)"
        )
    result = [list() for _ in range(num_clusters)]
    for cluster_num, clustered_elem in zip(cluster_assignment, value_pool):
        result[cluster_num].extend(clustered_elem)
    return result


@dataclass
class ClusteringResult(object):
    """
    Stores the result of clustering as a ClusteredIDs object.
    If clustering was not meaningful (e.g. one sequence per cluster),
    a set of sequences can be directly used as set of alternative
    alleles of the variant site under construction.
    """
    clustered_ids: ClusteredIDs
    sequences: Optional[Sequences] = None

    @property
    def no_clustering(self):
        return len(self.clustered_ids) == 1

    @property
    def have_precomputed_sequences(self):
        return self.sequences is not None

    def __str__(self):
        return self.__repr__()


def merge_sequences(*seqlists: Sequences, first_seq: str) -> Sequences:
    first_seq_found = False
    result = list()
    for seqlist in seqlists:
        for sequence in seqlist:
            if sequence == first_seq:
                first_seq_found = True
            else:
                result.append(sequence)

    # dev's responsibility to pass the correct first_seq argument
    assert first_seq_found, f"Provided first sequence argument ({first_seq}) not found in provided list of sequences " \
                            f"({seqlists})"

    merged_unexpanded_sequences = [first_seq] + result
    merged_expanded_sequences = SequenceExpander.get_expanded_sequences(merged_unexpanded_sequences)

    first_merged_expanded_sequence = merged_expanded_sequences[0]
    first_seq_has_ambiguous_char = first_seq != first_merged_expanded_sequence
    if first_seq_has_ambiguous_char:
        logger.warning(f"Provided first sequence argument ({first_seq}) has an ambiguous base, and thus does not match "
                       f"expanded first sequence ({first_merged_expanded_sequence})")

    return merged_expanded_sequences


def merge_clusters(*clusters: ClusteredIDs, first_id: str) -> ClusteredIDs:
    merged_clusters = list()
    first_id_cluster = []
    for cluster in chain.from_iterable(clusters):
        if first_id in cluster:
            first_id_cluster = cluster
        else:
            merged_clusters.append(cluster)
    if len(first_id_cluster) == 0:
        raise ValueError(f"Could not find {first_id} in any cluster")

    # place first_id in the first position of first_id_cluster
    first_id_cluster.remove(first_id)
    first_id_cluster.insert(0, first_id)
    return [first_id_cluster] + merged_clusters


def kmeans_cluster_seqs(
    alignment: MSA,
    kmer_size: int,
) -> ClusteringResult:
    """Divide sequences into subgroups of similar sequences.
    If no meaningful clustering is found, returns the (deduplicated, ungapped)
    set of input sequences.
    """
    # Find unique sequences for clustering, but keep each sequence's IDs
    seq_to_ids: SeqToIDs = defaultdict(list)
    seq_to_gapped_seqs: SeqToSeqs = defaultdict(list)
    small_seq_to_ids: SeqToIDs = defaultdict(list)
    first_id = alignment[0].id
    first_sequence = ungap(str(alignment[0].seq))

    for record in alignment:
        seq_with_gaps = str(record.seq)
        seq = ungap(seq_with_gaps)
        if len(seq) >= kmer_size:
            seq_to_ids[seq].append(record.id)
            seq_to_gapped_seqs[seq].append(seq_with_gaps)
        else:
            small_seq_to_ids[seq].append(record.id)

    num_clusters = 1
    num_sequences = len(seq_to_ids)
    too_few_seqs_to_cluster = num_sequences <= 2
    if too_few_seqs_to_cluster:
        single_cluster = [flatten_list(seq_to_ids.values()) + flatten_list(small_seq_to_ids.values())]
        clustered_ids = merge_clusters(single_cluster, first_id=first_id)
        merged_sequences = merge_sequences(seq_to_ids.keys(), small_seq_to_ids.keys(), first_seq=first_sequence)
        return ClusteringResult(clustered_ids, merged_sequences)

    distinct_sequences = list(seq_to_ids)
    distinct_kmers = count_distinct_kmers(distinct_sequences, kmer_size)
    count_matrix = count_kmer_occurrences(distinct_sequences, distinct_kmers)
    cluster_assignment = [0 for _ in range(len(seq_to_ids))]
    seqclustering: ClusteredSeqs = extract_clusters(
        seq_to_gapped_seqs, cluster_assignment
    )

    while cluster_further(seqclustering):
        num_clusters += 1
        if num_clusters > MAX_CLUSTERS:
            break
        if num_clusters == num_sequences:
            break
        kmeans = KMeans(n_clusters=num_clusters, random_state=2).fit(count_matrix)
        prev_cluster_assignment = cluster_assignment
        cluster_assignment = list(kmeans.predict(count_matrix))
        num_fitted_clusters = len(set(cluster_assignment))
        # Below holds when alignments are different, but kmer counts are identical
        # (due to repeats), making kmeans unable to fit requested number of clusters
        if num_fitted_clusters < num_clusters:
            cluster_assignment = prev_cluster_assignment
            num_clusters -= 1
            break
        seqclustering = extract_clusters(seq_to_gapped_seqs, cluster_assignment)

    no_clustering = num_clusters == 1 or num_clusters == num_sequences
    if no_clustering:
        single_cluster = [flatten_list(seq_to_ids.values()) + flatten_list(small_seq_to_ids.values())]
        clustered_ids = merge_clusters(single_cluster, first_id=first_id)
        merged_sequences = merge_sequences(seq_to_ids.keys(), small_seq_to_ids.keys(), first_seq=first_sequence)
        return ClusteringResult(clustered_ids, merged_sequences)
    else:
        clustered_ids: ClusteredIDs = []
        if num_sequences > 0:
            clustered_ids = extract_clusters(seq_to_ids, cluster_assignment)
        clustered_ids = merge_clusters(clustered_ids, small_seq_to_ids.values(), first_id=first_id)
        assert len(alignment) == sum(
            [len(i) for i in clustered_ids]
        ), "Each input sequence should be in a cluster"
        return ClusteringResult(clustered_ids)
