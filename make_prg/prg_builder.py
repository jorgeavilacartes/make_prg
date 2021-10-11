from typing import Tuple, Dict, Optional
from make_prg.io_utils import load_alignment_file
import pickle
from pathlib import Path
from make_prg.msa_aligner import MSAAligner
from make_prg.recursion_tree import SingleClusterNode, RecursiveTreeNode

class LeafNotFoundException(Exception):
    pass


class PrgBuilder(object):
    """
    Prg builder based from a multiple sequence alignment.
    """
    def __init__(
        self,
        locus_name: str,
        msa_file: Path,
        alignment_format: str,
        max_nesting: int,
        min_match_length: int,
        aligner: Optional[MSAAligner] = None,
    ):
        self.locus_name: str = locus_name
        self.max_nesting: int = max_nesting
        self.min_match_length: int = min_match_length
        self.aligner: Optional[MSAAligner] = aligner
        self.next_node_id: int = 0
        self._site_num: int = 5
        self.leaves_index: Dict[Tuple[int, int], RecursiveTreeNode] = {}

        alignment = load_alignment_file(str(msa_file), alignment_format)
        self.root: RecursiveTreeNode = SingleClusterNode(
            nesting_level=0,
            alignment=alignment,
            parent=None,
            prg_builder=self
        )

    def build_prg(self) -> str:
        self._site_num = 5
        prg_as_list = []
        self.root.preorder_traversal_to_build_prg(prg_as_list)
        prg = "".join(prg_as_list)
        return prg

    def get_next_site_num(self) -> int:
        previous_site_num = self._site_num
        self._site_num += 2
        return previous_site_num

    def get_next_node_id(self) -> int:
        self.next_node_id += 1
        return self.next_node_id - 1

    def update_leaves_index(self, start_index: int, end_index: int, node: RecursiveTreeNode):
        interval = (start_index, end_index)
        self.leaves_index[interval] = node

    def get_node_given_interval(self, interval: Tuple[int, int]) -> RecursiveTreeNode:
        # TODO: move this back to assert once is solved
        interval_is_indexed = interval in self.leaves_index
        if not interval_is_indexed:
            raise LeafNotFoundException(
                f"Queried interval {interval} does not exist in leaves index for locus {self.locus_name}"
            )

        # assert interval in self.leaves_index, \
        #     f"Fatal error: Queried interval {interval} does not exist in leaves index for locus {self.locus_name}"

        return self.leaves_index[interval]

    def serialize(self, filepath: [Path, str]):
        with open(filepath, "wb") as filehandler:
            pickle.dump(self, filehandler)

    @staticmethod
    def deserialize(filepath: [Path, str]) -> "PrgBuilder":
        with open(filepath, "rb") as filehandler:
            return pickle.load(filehandler)


class PrgBuilderCollection:
    """
    Represent a collection of PrgBuilder, to be serialised and deserialised
    """
    def __init__(self, locus_name_to_pickle_files: Dict[str, str]):
        self.locus_name_to_pickle_files: Dict[str, str] = locus_name_to_pickle_files

    def serialize(self, filepath: [Path, str]):
        with open(filepath, "wb") as filehandler:
            pickle.dump(self, filehandler)

    @staticmethod
    def deserialize(filepath: [Path, str]) -> "PrgBuilderCollection":
        with open(filepath, "rb") as filehandler:
            return pickle.load(filehandler)

    def to_absolute_paths_wrt_given_parent(self, parent: Path):
        for locus_name, pickle_file in self.locus_name_to_pickle_files.items():
            absolute_filepath = parent / pickle_file
            self.locus_name_to_pickle_files[locus_name] = str(absolute_filepath)
