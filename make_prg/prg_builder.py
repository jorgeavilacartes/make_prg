from typing import Tuple

from make_prg.io_utils import load_alignment_file
import pickle
from pathlib import Path
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from make_prg.msa_aligner import MSAAligner
from make_prg.recursion_tree import SingleClusterNode


class LeafNotFoundException(Exception):
    pass


class PrgBuilder(object):
    """
    Prg builder based from a multiple sequence alignment.
    Note min_match_length must be strictly greater than max_nesting + 1. TODO: why?
    """
    def __init__(
        self,
        locus_name,
        msa_file,
        alignment_format,
        max_nesting,
        min_match_length,
        aligner: MSAAligner,
    ):
        self.locus_name = locus_name
        self.max_nesting = max_nesting
        self.min_match_length = min_match_length
        self.aligner = aligner
        self.next_node_id = 0
        self.leaves_index = {}

        alignment = load_alignment_file(msa_file, alignment_format)
        self._root = SingleClusterNode(
            nesting_level=1,
            alignment=alignment,
            parent=None,
            prg_builder=self
        )

    def build_prg(self):
        self._site_num = 5
        prg_as_list = []
        self._root.preorder_traversal_to_build_prg(prg_as_list)
        prg = "".join(prg_as_list)
        return prg

    def get_next_site_num(self):
        previous_site_num = self._site_num
        self._site_num += 2
        return previous_site_num

    def get_next_node_id(self):
        self.next_node_id += 1
        return self.next_node_id - 1

    def update_leaves_index(self, start_index: int, end_index: int, node):
        interval = (start_index, end_index)
        self.leaves_index[interval] = node

    def get_node_given_interval(self, interval: Tuple[int, int]):
        # TODO: move this back to assert once is solved
        interval_is_indexed = interval in self.leaves_index
        if not interval_is_indexed:
            raise LeafNotFoundException(
                f"Queried interval {interval} does not exist in leaves index for locus {self.locus_name}"
            )

        # assert interval in self.leaves_index, \
        #     f"Fatal error: Queried interval {interval} does not exist in leaves index for locus {self.locus_name}"

        return self.leaves_index[interval]

    def serialize(self, filename):
        with open(filename, "wb") as filehandler:
            pickle.dump(self, filehandler)

    @staticmethod
    def deserialize(filename):
        with open(filename, "rb") as filehandler:
            return pickle.load(filehandler)

    def output_graph(self, filename):
        drawing_level_to_node_ids = defaultdict(list)
        for node in self.all_nodes:
            drawing_level_to_node_ids[node.drawing_level].append(node.id)

        plt.figure(figsize=(20, 10))
        a_graph = nx.drawing.nx_agraph.to_agraph(self.recursion_tree)
        for node in a_graph.nodes():
            node.attr["label"] = self.leaf_id_to_seq.get(int(node.name), node.name)

        for drawing_level, nodes in drawing_level_to_node_ids.items():
            a_graph.add_subgraph(nodes, rank="same")
        a_graph.layout(prog="dot")
        a_graph.draw(filename)
        a_graph.draw(f"{filename}.dot")


class PrgBuilderCollection:
    """
    Represent a collection of PrgBuilder, to be serialised and deserialised
    """
    def __init__(self, locus_name_to_pickle_files):
        self.locus_name_to_pickle_files = locus_name_to_pickle_files

    def serialize(self, filepath: [Path, str]):
        with open(filepath, "wb") as filehandler:
            pickle.dump(self, filehandler)

    @staticmethod
    def deserialize(filepath: [Path, str]):
        with open(filepath, "rb") as filehandler:
            return pickle.load(filehandler)

    def to_absolute_paths_wrt_given_parent(self, parent: Path):
        for locus_name, pickle_file in self.locus_name_to_pickle_files.items():
            absolute_filepath = parent / pickle_file
            self.locus_name_to_pickle_files[locus_name] = str(absolute_filepath)
