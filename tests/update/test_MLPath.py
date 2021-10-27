from unittest import TestCase
from make_prg.update.MLPath import MLPathNode, MLPath
from intervaltree import IntervalTree, Interval

class MLPathTest(TestCase):
    def setUp(self):
        self.ml_path_node_1 = MLPathNode(key=(5, 9), sequence="ACGT")
        self.ml_path_node_2 = MLPathNode(key=(12, 13), sequence="G")
        self.ml_path_node_3 = MLPathNode(key=(20, 22), sequence="AA")

    # TODO: might not be the best practice to check the private attributes in tests
    # TODO: it breaks encapsulation
    def test___constructor___single_MLPath(self):
        ml_path_nodes = [self.ml_path_node_1]
        ml_path = MLPath(ml_path_nodes)

        self.assertEqual(ml_path_nodes, ml_path._ml_path_nodes)

        expected_index = IntervalTree([Interval(0, 4, MLPathNode(key=(5, 9), sequence="ACGT"))])
        self.assertEqual(expected_index, ml_path._ml_path_index)

    def test___constructor___MLPath_with_two_nodes(self):
        ml_path_nodes = [self.ml_path_node_1, self.ml_path_node_2]
        ml_path = MLPath(ml_path_nodes)

        self.assertEqual(ml_path_nodes, ml_path._ml_path_nodes)

        expected_index = IntervalTree([
            Interval(0, 4, MLPathNode(key=(5, 9), sequence="ACGT")),
            Interval(4, 5, MLPathNode(key=(12, 13), sequence="G"))
        ])
        self.assertEqual(expected_index, ml_path._ml_path_index)

    def test___constructor___MLPath_with_three_nodes(self):
        ml_path_nodes = [self.ml_path_node_1, self.ml_path_node_2, self.ml_path_node_3]
        ml_path = MLPath(ml_path_nodes)

        self.assertEqual(ml_path_nodes, ml_path._ml_path_nodes)

        expected_index = IntervalTree([
            Interval(0, 4, MLPathNode(key=(5, 9), sequence="ACGT")),
            Interval(4, 5, MLPathNode(key=(12, 13), sequence="G")),
            Interval(5, 7, MLPathNode(key=(20, 22), sequence="AA"))
        ])
        self.assertEqual(expected_index, ml_path._ml_path_index)
