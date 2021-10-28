from unittest import TestCase
from make_prg.update.MLPath import MLPathNode, MLPath, MLPathError
from intervaltree import IntervalTree, Interval


class MLPathTest(TestCase):
    def setUp(self):
        self.ml_path_node_1 = MLPathNode(key=(5, 9), sequence="ACGT")
        self.ml_path_node_2 = MLPathNode(key=(12, 13), sequence="G")
        self.ml_path_node_3 = MLPathNode(key=(20, 22), sequence="AA")
        self.ml_path_nodes = [self.ml_path_node_1, self.ml_path_node_2, self.ml_path_node_3]
        self.ml_path = MLPath(self.ml_path_nodes)

    def test___constructor___empty_MLPath___raises_MLPathError(self):
        with self.assertRaises(MLPathError):
            MLPath([])

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
        self.assertEqual(self.ml_path_nodes, self.ml_path._ml_path_nodes)
        expected_index = IntervalTree([
            Interval(0, 4, MLPathNode(key=(5, 9), sequence="ACGT")),
            Interval(4, 5, MLPathNode(key=(12, 13), sequence="G")),
            Interval(5, 7, MLPathNode(key=(20, 22), sequence="AA"))
        ])
        self.assertEqual(expected_index, self.ml_path._ml_path_index)

    def test___get_node_at_position___positions_of_ml_path_node_1(self):
        expected = self.ml_path_node_1
        for position in range(0, 4):
            actual = self.ml_path.get_node_at_position(position)
            self.assertEqual(expected, actual)

    def test___get_node_at_position___positions_of_ml_path_node_2(self):
        expected = self.ml_path_node_2
        actual = self.ml_path.get_node_at_position(4)
        self.assertEqual(expected, actual)

    def test___get_node_at_position___positions_of_ml_path_node_3(self):
        expected = self.ml_path_node_3
        for position in range(5, 7):
            actual = self.ml_path.get_node_at_position(position)
            self.assertEqual(expected, actual)

    def test___get_node_at_position___one_position_over_the_last_one___last_insert_position(self):
        expected = self.ml_path_node_3
        actual = self.ml_path.get_node_at_position(7)
        self.assertEqual(expected, actual)

    def test___get_node_at_position___two_positions_over_the_last_one___raises_MLPathError(self):
        with self.assertRaises(MLPathError):
            self.ml_path.get_node_at_position(8)
