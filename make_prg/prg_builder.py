from typing import Tuple, Dict, Optional, List
from make_prg.utils.io_utils import load_alignment_file, SetOutputFiles
import pickle
from pathlib import Path
from zipfile import ZipFile
from make_prg.utils.msa_aligner import MSAAligner
from make_prg.recursion_tree import SingleClusterNode, RecursiveTreeNode
from make_prg.utils.drawer import RecursiveTreeDrawer
import os


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
        self.site_num: int = 5
        self.leaves_index: Dict[Tuple[int, int], SingleClusterNode] = {}

        alignment = load_alignment_file(str(msa_file), alignment_format)
        self.root: RecursiveTreeNode = SingleClusterNode(
            nesting_level=0,
            alignment=alignment,
            parent=None,
            prg_builder=self
        )

    def build_prg(self) -> str:
        self.site_num = 5
        prg_as_list = []
        self.root.preorder_traversal_to_build_prg(prg_as_list)
        prg = "".join(prg_as_list)
        return prg

    def get_next_site_num(self) -> int:
        site_num = self.site_num
        self.site_num += 2
        return site_num

    def get_next_node_id(self) -> int:
        self.next_node_id += 1
        return self.next_node_id - 1

    def update_leaves_index(self, start_index: int, end_index: int, node: SingleClusterNode):
        interval = (start_index, end_index)
        self.leaves_index[interval] = node

    def get_node_given_interval(self, interval: Tuple[int, int]) -> SingleClusterNode:
        # TODO: move this back to assert once is solved
        # TODO: should it really be an assert?
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
    def deserialize(bytes_from_zip) -> "PrgBuilder":
        return pickle.loads(bytes_from_zip)

    def output_debug_graphs(self, debug_graphs_dir):
        os.makedirs(debug_graphs_dir, exist_ok=True)
        recursive_tree_drawer = RecursiveTreeDrawer(self.root)
        recursive_tree_drawer.output_graph(debug_graphs_dir / f"{self.locus_name}.recursion_tree.png")


class PrgBuilderZipDatabase:
    """
    Represent a collection of PrgBuilders, to be saved to and loaded from a zip file
    """
    def __init__(self, zip_filepath: Path):
        is_a_zip_file = zip_filepath.suffix == ".zip"
        assert is_a_zip_file, "PrgBuilderZipDatabase initialised without a .zip filepath"
        self._zip_filepath: Path = zip_filepath
        self._zip_file: Optional[ZipFile] = None

    def save(self, locus_to_prg_builder_pickle_path: Dict[str, Path]):
        with ZipFile(self._zip_filepath, "w") as zip_file:
            for locus, prg_builder_pickle_path in locus_to_prg_builder_pickle_path.items():
                zip_file.write(prg_builder_pickle_path, locus)

    def load(self):
        self._zip_file = ZipFile(self._zip_filepath)

    def close(self):
        if self._zip_file is not None:
            self._zip_file.close()

    def get_loci_names(self) -> List[str]:
        return self._zip_file.namelist()

    def get_PrgBuilder(self, locus: str) -> PrgBuilder:
        return PrgBuilder.deserialize(self._zip_file.read(locus))
