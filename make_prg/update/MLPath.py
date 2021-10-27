from typing import Optional, List, Tuple
from intervaltree.intervaltree import IntervalTree


class MLPathError(Exception):
    pass


class EmptyMLPathSequence(Exception):
    pass


class MLPathNode:
    def __init__(self, key: Tuple[int, int], sequence: str):
        self.key: Tuple[int, int] = key
        self._set_sequence(sequence)

        # these are set by class MLPath during indexing
        self.start_index_in_linear_path: Optional[int] = None
        self.end_index_in_linear_path: Optional[int] = None

        self._check_is_a_valid_node()

    def _set_sequence(self, sequence: str):
        empty_ML_path_sequence = len(sequence) == 0
        if empty_ML_path_sequence:
            raise EmptyMLPathSequence(f"Found a ML path node ({self.key}) with empty sequence")
        self.sequence: str = sequence

    def _check_is_a_valid_node(self):
        interval_size = self.key[1] - self.key[0]
        sequence_size = len(self.sequence)
        valid_node = interval_size == sequence_size
        if not valid_node:
            raise MLPathError(f"{self} is not a valid node")

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return (self.key, self.sequence) == (other.key, other.sequence)
        else:
            return False

    def __str__(self):
        return f"PRG key = {self.key}; " \
               f"ML seq interval = [{self.start_index_in_linear_path}:{self.end_index_in_linear_path}]; " \
               f"Seq = {self.sequence}"

    def __repr__(self):
        return f"MLPathNode(key={self.key}, sequence=\"{self.sequence}\")"


class MLPath:
    def __init__(self, ml_path_nodes: List[MLPathNode]):
        if len(ml_path_nodes) == 0:
            raise MLPathError("ML paths cannot be empty")
        self._ml_path_nodes: List[MLPathNode] = ml_path_nodes
        self._ml_path_index: IntervalTree = IntervalTree()
        self._index()

    def _index(self):
        start_index_in_linear_path = 0
        for ml_path_node_index, ml_path_node in enumerate(self._ml_path_nodes):
            end_index_in_linear_path = start_index_in_linear_path + len(ml_path_node.sequence)
            self._ml_path_index.addi(
                start_index_in_linear_path,
                end_index_in_linear_path,
                data=ml_path_node,
            )
            ml_path_node.start_index_in_linear_path = start_index_in_linear_path
            ml_path_node.end_index_in_linear_path = end_index_in_linear_path
            start_index_in_linear_path = end_index_in_linear_path

    def _last_insertion_pos(self):
        """
        This position is not indexed, but can be one where we can insert a sequence in the last node
        """
        return self._ml_path_nodes[-1].end_index_in_linear_path

    def get_node_at_position(self, position: int) -> MLPathNode:
        is_an_insertion_in_the_last_position_of_the_last_node = position == self._last_insertion_pos()
        if is_an_insertion_in_the_last_position_of_the_last_node:
            last_node = self._ml_path_nodes[-1]
            return last_node

        # now to the common case
        nodes = self._ml_path_index[position]

        two_or_more_nodes_overlap_this_position = len(nodes) >= 2
        assert not two_or_more_nodes_overlap_this_position, f"2+ nodes overlap at position {position}, " \
                                                            f"this ML path ({self}) probably has an indexing bug."
        no_nodes_overlap_this_position = len(nodes) == 0
        if not no_nodes_overlap_this_position:
            # this is an assertion, as this can be an input error (e.g. a badly formed denovo_paths.txt file)
            raise MLPathError(f"No nodes overlap this ML path ({self}) at position {position}, "
                              f"is the denovo_paths.txt given as input correct?")

        # only one node overlap this position
        node = list(nodes)[0].data
        return node
