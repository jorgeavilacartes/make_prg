from unittest import TestCase
from unittest.mock import patch, Mock, PropertyMock
from tests.test_helpers import make_alignment, equal_msas, first_dict_contained_in_second
from make_prg.recursion_tree import MultiClusterNode, SingleClusterNode
from make_prg.prg_builder import PrgBuilder
from make_prg.from_msa import MSA
from pathlib import Path
from make_prg.from_msa.interval_partition import IntervalType, Interval


@patch.object(PrgBuilder, PrgBuilder.get_next_node_id.__name__, return_value=0)
@patch.object(SingleClusterNode, PrgBuilder.__init__.__name__, return_value=None)
@patch("make_prg.prg_builder.load_alignment_file")
class TestMultiClusterNode(TestCase):
    # Note: can't apply patches to setUp(), so creating this method that is called in every test
    def setup(self) -> None:
        self.alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG"], ["s1", "s2", "s3", "s4"]
        )
        self.prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        self.multi_cluster_node = MultiClusterNode(1, self.alignment, None, self.prg_builder, False)

    @patch.object(MultiClusterNode, MultiClusterNode._get_children.__name__, return_value=["child_1", "child_2"])
    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA", return_value="remove_columns_full_of_gaps_from_MSA_mock")
    def test___constructor(self, *uninteresting_mocks):
        self.setup()
        parent_mock = Mock()
        node = MultiClusterNode(1, self.alignment, parent_mock, self.prg_builder, False)

        self.assertEqual(1, node.nesting_level)
        self.assertTrue("remove_columns_full_of_gaps_from_MSA_mock", node.alignment)
        self.assertEqual(parent_mock, node.parent)
        self.assertEqual(self.prg_builder, node.prg_builder)
        self.assertFalse(node.force_no_child)
        self.assertEqual(0, node.id)
        self.assertEqual(set(), node.new_sequences)
        self.assertEqual(["child_1", "child_2"], node.children)

    @patch.object(MultiClusterNode, MultiClusterNode._get_children.__name__, return_value=[])
    def test___constructor___no_children___raises_AssertionError(self, *uninteresting_mocks):
        alignment = make_alignment(["A"])
        prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        with self.assertRaises(AssertionError):
            MultiClusterNode(1, alignment, None, prg_builder, False)




    class Preorder_traversal_to_build_prg_Mock:
        def __init__(self, child_id):
            self.child_id = child_id
        def preorder_traversal_to_build_prg(self, prg_as_list, delim_char):
            prg_as_list.append(f"child_{self.child_id}")

    @patch.object(PrgBuilder, PrgBuilder.get_next_site_num.__name__, return_value=32)
    def test___preorder_traversal_to_build_prg___single_child(self, *mocks):
        self.setup()
        node = MultiClusterNode(1, self.alignment, None, self.prg_builder, False)

        # get the original method
        preorder_traversal_to_build_prg = MultiClusterNode.preorder_traversal_to_build_prg

        # now patch the children
        children = [TestMultiClusterNode.Preorder_traversal_to_build_prg_Mock(1)]
        with patch.object(MultiClusterNode, "children", new_callable=PropertyMock, return_value=children):
            # call the original method with the recursive calls now patched
            prg_as_list = []
            preorder_traversal_to_build_prg(node, prg_as_list, delim_char="*")

            expected_prg = "*32*child_1*32*"
            actual_prg = "".join(prg_as_list)

            self.assertEqual(expected_prg, actual_prg)

    @patch.object(PrgBuilder, PrgBuilder.get_next_site_num.__name__, return_value=32)
    def test___preorder_traversal_to_build_prg___multiple_children(self, *mocks):
        self.setup()
        node = MultiClusterNode(1, self.alignment, None, self.prg_builder, False)

        # get the original method
        preorder_traversal_to_build_prg = MultiClusterNode.preorder_traversal_to_build_prg

        # now patch the children
        children = [TestMultiClusterNode.Preorder_traversal_to_build_prg_Mock(1),
                    TestMultiClusterNode.Preorder_traversal_to_build_prg_Mock(2),
                    TestMultiClusterNode.Preorder_traversal_to_build_prg_Mock(3)]
        with patch.object(MultiClusterNode, "children", new_callable=PropertyMock, return_value=children):
            # call the original method with the recursive calls now patched
            prg_as_list = []
            preorder_traversal_to_build_prg(node, prg_as_list, delim_char="*")

            expected_prg = "*32*child_1*33*child_2*33*child_3*32*"
            actual_prg = "".join(prg_as_list)

            self.assertEqual(expected_prg, actual_prg)



    @patch("make_prg.recursion_tree.kmeans_cluster_seqs", return_value=[["s5", "s0", "s1", "s6"], ["s3"], ["s2", "s4"]])
    def test___get_subalignments_by_clustering(self, *uninteresting_mocks):
        alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG", "CCCC", "TTTT", "AAAA"], ["s0", "s1", "s2", "s3", "s4", "s5", "s6"]
        )
        prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        node = MultiClusterNode(1, alignment, None, prg_builder, False)

        expected = [
            make_alignment(
                ["AAAT", "C--C", "TTTT", "AAAA"], ["s0", "s1", "s5", "s6"]),
            make_alignment(
                ["GNGG"], ["s3"]),
            make_alignment(
                ["AATT", "CCCC"],
                ["s2", "s4"]),
        ]
        actual = node._get_subalignments_by_clustering()

        self.assertEqual(3, len(actual))
        for i in range(3):
            self.assertTrue(equal_msas(expected[i], actual[i]))

    def test___get_sub_alignment_by_list_id___GivenOrderedIds_SubalignmentInSequenceOrder(self, *uninteresting_mocks):
        self.setup()
        expected = MSA([self.alignment[0], self.alignment[2]])
        actual = self.multi_cluster_node._get_sub_alignment_by_list_id(["s1", "s3"])
        self.assertTrue(equal_msas(expected, actual))

    def test___get_sub_alignment_by_list_id___GivenUnorderedIds_SubalignmentStillInSequenceOrder(self, *uninteresting_mocks):
        """
        Sequences given rearranged are still output in input order
        """
        self.setup()
        expected = MSA([self.alignment[0], self.alignment[2]])
        actual = self.multi_cluster_node._get_sub_alignment_by_list_id(["s3", "s1"])
        self.assertTrue(equal_msas(expected, actual))

    @patch("make_prg.recursion_tree.SingleClusterNode")
    @patch.object(MultiClusterNode, MultiClusterNode._get_subalignments_by_clustering.__name__)
    def test___get_children___single_cluster(self, get_subalignments_by_clustering_mock, SingleClusterNode_mock,
                                             *mocks):
        single_alignment = make_alignment(["AAAT", "C--C", "AATT"])
        subalignments = [single_alignment]
        get_subalignments_by_clustering_mock.return_value = subalignments
        self.setup()

        expected = [SingleClusterNode_mock.return_value]
        actual = self.multi_cluster_node.children

        self.assertEqual(expected, actual)
        SingleClusterNode_mock.assert_called_once_with(
            nesting_level=1,
            alignment=single_alignment,
            parent=self.multi_cluster_node,
            prg_builder=self.prg_builder,
            force_no_child=True
        )

    get_children___multiple_clusters___SingleClusterNode___side_effect = [Mock(), Mock(), Mock()]
    @patch("make_prg.recursion_tree.SingleClusterNode", side_effect=get_children___multiple_clusters___SingleClusterNode___side_effect)
    @patch.object(MultiClusterNode, MultiClusterNode._get_subalignments_by_clustering.__name__)
    def test___get_children___multiple_clusters(self, get_subalignments_by_clustering_mock, SingleClusterNode_mock,
                                             *mocks):
        alignment_1 = make_alignment(["AAAT", "C--C", "AATT"])
        alignment_2 = make_alignment(["GGGG"])
        alignment_3 = make_alignment(["CCCC", "TTTT"])

        subalignments = [alignment_1, alignment_2, alignment_3]
        get_subalignments_by_clustering_mock.return_value = subalignments
        self.setup()

        expected = self.get_children___multiple_clusters___SingleClusterNode___side_effect
        actual = self.multi_cluster_node.children

        self.assertEqual(expected, actual)
        self.assertEqual(3, SingleClusterNode_mock.call_count)
        SingleClusterNode_mock.assert_any_call(
            nesting_level=1,
            alignment=alignment_1,
            parent=self.multi_cluster_node,
            prg_builder=self.prg_builder,
            force_no_child=False
        )
        SingleClusterNode_mock.assert_any_call(
            nesting_level=1,
            alignment=alignment_2,
            parent=self.multi_cluster_node,
            prg_builder=self.prg_builder,
            force_no_child=False
        )
        SingleClusterNode_mock.assert_any_call(
            nesting_level=1,
            alignment=alignment_3,
            parent=self.multi_cluster_node,
            prg_builder=self.prg_builder,
            force_no_child=False
        )


@patch.object(PrgBuilder, PrgBuilder.get_next_node_id.__name__, return_value=0)
@patch("make_prg.prg_builder.load_alignment_file")
class TestSingleClusterNode(TestCase):
    def subsetup(self):
        self.alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG"], ["s1", "s2", "s3", "s4"]
        )
        with patch.object(SingleClusterNode, PrgBuilder.__init__.__name__, return_value=None):
            self.prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)

    # Note: can't apply patches to setUp(), so creating this method that is called in every test
    def setup(self) -> None:
        self.subsetup()
        self.single_cluster_node = SingleClusterNode(1, self.alignment, None, self.prg_builder, False)

    @patch.object(SingleClusterNode, SingleClusterNode._get_children.__name__, return_value=["child_1", "child_2"])
    @patch("make_prg.recursion_tree.remove_columns_full_of_gaps_from_MSA",
           return_value="remove_columns_full_of_gaps_from_MSA_mock")
    @patch("make_prg.recursion_tree.get_consensus_from_MSA", return_value="consensus_from_MSA")
    @patch("make_prg.recursion_tree.IntervalPartitioner")
    def test___constructor(self, IntervalPartitioner_mock, get_consensus_from_MSA_mock, *uninteresting_mocks):
        alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG"], ["s1", "s2", "s3", "s4"]
        )
        with patch.object(SingleClusterNode, PrgBuilder.__init__.__name__, return_value=None):
            prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        parent_mock = Mock()
        class GetIntervalsMock:
            def get_intervals(self):
                return "match_intervals", "non_match_intervals", "all_intervals"
        get_intervals_mock = GetIntervalsMock()
        IntervalPartitioner_mock.return_value = get_intervals_mock

        node = SingleClusterNode(1, alignment, parent_mock, prg_builder, False)

        self.assertEqual(1, node.nesting_level)
        self.assertTrue("remove_columns_full_of_gaps_from_MSA_mock", node.alignment)
        self.assertEqual(parent_mock, node.parent)
        self.assertEqual(prg_builder, node.prg_builder)
        self.assertFalse(node.force_no_child)
        self.assertEqual(0, node.id)
        self.assertEqual(set(), node.new_sequences)
        self.assertEqual(["child_1", "child_2"], node.children)
        self.assertEqual("consensus_from_MSA", node.consensus)
        self.assertEqual(len("consensus_from_MSA"), node.length)
        self.assertEqual("match_intervals", node.match_intervals)
        self.assertEqual("non_match_intervals", node.non_match_intervals)
        self.assertEqual("all_intervals", node.all_intervals)
        get_consensus_from_MSA_mock.assert_called_once_with("remove_columns_full_of_gaps_from_MSA_mock")
        IntervalPartitioner_mock.assert_called_once_with("consensus_from_MSA", 7, "remove_columns_full_of_gaps_from_MSA_mock")

    @patch.object(SingleClusterNode, SingleClusterNode._get_children.__name__, return_value=[])
    def test___is_leaf(self, *uninteresting_mocks):
        self.setup()
        node = SingleClusterNode(1, self.alignment, None, self.prg_builder, False)

        self.assertTrue(node.is_leaf())

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=1)
    def test___alignment_has_issues___num_unique_nongapped_seqs_1(self,
                 get_number_of_unique_ungapped_sequences_mock, *uninteresting_mocks):
        self.setup()
        self.assertTrue(self.single_cluster_node._alignment_has_issues(self.single_cluster_node.alignment))
        get_number_of_unique_ungapped_sequences_mock.assert_called_once_with(self.single_cluster_node.alignment)

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=2)
    def test___alignment_has_issues___num_unique_nongapped_seqs_2(self,
                 get_number_of_unique_ungapped_sequences_mock, *uninteresting_mocks):
        self.setup()
        self.assertTrue(self.single_cluster_node._alignment_has_issues(self.single_cluster_node.alignment))
        get_number_of_unique_ungapped_sequences_mock.assert_called_once_with(self.single_cluster_node.alignment)

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=3)
    @patch("make_prg.recursion_tree.get_number_of_unique_gapped_sequences", return_value=4)
    def test___alignment_has_issues___alignment_has_ambiguity___ungapped_smaller_gapped(self,
         get_number_of_unique_gapped_sequences_mock, get_number_of_unique_ungapped_sequences_mock, *uninteresting_mocks):
        self.setup()
        self.assertTrue(self.single_cluster_node._alignment_has_issues(self.single_cluster_node.alignment))
        get_number_of_unique_gapped_sequences_mock.assert_called_once_with(self.single_cluster_node.alignment)
        get_number_of_unique_ungapped_sequences_mock.assert_called_once_with(self.single_cluster_node.alignment)

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=4)
    @patch("make_prg.recursion_tree.get_number_of_unique_gapped_sequences", return_value=4)
    def test___alignment_has_issues___alignment_has_ambiguity___ungapped_equals_gapped___no_issues(self,
        get_number_of_unique_gapped_sequences_mock, get_number_of_unique_ungapped_sequences_mock, *uninteresting_mocks):
        self.setup()
        self.assertFalse(self.single_cluster_node._alignment_has_issues(self.single_cluster_node.alignment))
        get_number_of_unique_gapped_sequences_mock.assert_called_once_with(self.single_cluster_node.alignment)
        get_number_of_unique_ungapped_sequences_mock.assert_called_once_with(self.single_cluster_node.alignment)

    @patch("make_prg.recursion_tree.get_number_of_unique_ungapped_sequences", return_value=5)
    @patch("make_prg.recursion_tree.get_number_of_unique_gapped_sequences", return_value=4)
    def test___alignment_has_issues___alignment_has_ambiguity___ungapped_larger_gapped___raises_AssertionError(self,
        *uninteresting_mocks):
        self.setup()
        with self.assertRaises(AssertionError):
            self.single_cluster_node._alignment_has_issues(self.single_cluster_node.alignment)

    def test___infer_if_this_node_should_have_no_child___force_no_child(self, *uninteresting_mocks):
        self.setup()
        self.single_cluster_node.force_no_child = True

        self.assertTrue(self.single_cluster_node._infer_if_this_node_should_have_no_child())

    def test___infer_if_this_node_should_have_no_child___single_match_interval(self, *uninteresting_mocks):
        self.setup()
        interval = Interval(IntervalType.Match, 3, 10)
        self.single_cluster_node.force_no_child = False
        self.single_cluster_node.all_intervals = [interval]
        self.single_cluster_node.match_intervals = [interval]

        self.assertTrue(self.single_cluster_node._infer_if_this_node_should_have_no_child())

    def test___infer_if_this_node_should_have_no_child___max_nesting_level_reached(self, *uninteresting_mocks):
        self.setup()
        self.single_cluster_node.force_no_child = False
        self.single_cluster_node.all_intervals = [Interval(IntervalType.Match, 0, 5),
                                                  Interval(IntervalType.Match, 6, 10)]
        self.single_cluster_node.nesting_level = 5

        self.assertTrue(self.single_cluster_node._infer_if_this_node_should_have_no_child())

    def test___infer_if_this_node_should_have_no_child___small_variant_site(self, *uninteresting_mocks):
        self.setup()
        self.single_cluster_node.force_no_child = False
        self.single_cluster_node.all_intervals = [Interval(IntervalType.Match, 0, 5),
                                                  Interval(IntervalType.Match, 6, 10)]
        self.single_cluster_node.nesting_level = 1
        self.single_cluster_node.prg_builder.min_match_length = 10

        self.assertTrue(self.single_cluster_node._infer_if_this_node_should_have_no_child())

    @patch.object(SingleClusterNode, SingleClusterNode._alignment_has_issues.__name__, return_value=True)
    def test___infer_if_this_node_should_have_no_child___alignment_has_issues(self,
            alignment_has_issues_mock, *uninteresting_mocks):
        self.setup()
        self.single_cluster_node.force_no_child = False
        self.single_cluster_node.all_intervals = [Interval(IntervalType.Match, 0, 5),
                                                  Interval(IntervalType.Match, 6, 10)]
        self.single_cluster_node.nesting_level = 1
        self.single_cluster_node.prg_builder.min_match_length = 1

        self.assertTrue(self.single_cluster_node._infer_if_this_node_should_have_no_child())
        alignment_has_issues_mock.assert_called_once_with(self.single_cluster_node.alignment)

    @patch.object(SingleClusterNode, SingleClusterNode._alignment_has_issues.__name__, return_value=False)
    def test___infer_if_this_node_should_have_no_child___no_issues___can_have_child(self,
            alignment_has_issues_mock, *uninteresting_mocks):
        self.setup()
        self.single_cluster_node.force_no_child = False
        self.single_cluster_node.all_intervals = [Interval(IntervalType.Match, 0, 5),
                                                  Interval(IntervalType.Match, 6, 10)]
        self.single_cluster_node.nesting_level = 1
        self.single_cluster_node.prg_builder.min_match_length = 1

        self.assertFalse(self.single_cluster_node._infer_if_this_node_should_have_no_child())
        alignment_has_issues_mock.assert_called_once_with(self.single_cluster_node.alignment)

    def test___infer_if_should_not_cluster___is_a_match_interval(self, *uninteresting_mocks):
        self.setup()
        interval = Interval(IntervalType.Match, 5, 10)
        self.single_cluster_node.match_intervals = [Interval(IntervalType.Match, 0, 5),
                                                    interval,
                                                    Interval(IntervalType.Match, 20, 25)]

        self.assertTrue(self.single_cluster_node._infer_if_should_not_cluster(interval, self.alignment))

    @patch.object(SingleClusterNode, SingleClusterNode._alignment_has_issues.__name__, return_value=True)
    def test___infer_if_should_not_cluster___alignment_has_issues(self,
          alignment_has_issues_mock, *uninteresting_mocks):
        self.setup()
        interval = Interval(IntervalType.Match, 5, 10)
        self.single_cluster_node.match_intervals = []

        self.assertTrue(self.single_cluster_node._infer_if_should_not_cluster(interval, self.single_cluster_node.alignment))
        alignment_has_issues_mock.assert_called_once_with(self.single_cluster_node.alignment)

    @patch("make_prg.recursion_tree.kmeans_cluster_seqs", return_value=[["s1", "s2", "s3", "s4"]])
    @patch.object(SingleClusterNode, SingleClusterNode._alignment_has_issues.__name__, return_value=False)
    def test___infer_if_should_not_cluster___single_cluster(self,
           alignment_has_issues_mock, kmeans_cluster_seqs_mock, *uninteresting_mocks):
        self.setup()
        interval = Interval(IntervalType.Match, 5, 10)
        self.single_cluster_node.match_intervals = []

        self.assertTrue(
            self.single_cluster_node._infer_if_should_not_cluster(interval, self.single_cluster_node.alignment))
        alignment_has_issues_mock.assert_called_once_with(self.single_cluster_node.alignment)
        kmeans_cluster_seqs_mock.assert_called_once_with(self.single_cluster_node.alignment, self.prg_builder.min_match_length)

    @patch("make_prg.recursion_tree.kmeans_cluster_seqs", return_value=[["s1", "s2"], ["s3", "s4"]])
    @patch.object(SingleClusterNode, SingleClusterNode._alignment_has_issues.__name__, return_value=False)
    def test___infer_if_should_not_cluster___two_clusters___no_issues_can_cluster(self,
                                                            alignment_has_issues_mock, kmeans_cluster_seqs_mock,
                                                            *uninteresting_mocks):
        self.setup()
        interval = Interval(IntervalType.Match, 5, 10)
        self.single_cluster_node.match_intervals = []

        self.assertFalse(
            self.single_cluster_node._infer_if_should_not_cluster(interval, self.single_cluster_node.alignment))
        alignment_has_issues_mock.assert_called_once_with(self.single_cluster_node.alignment)
        kmeans_cluster_seqs_mock.assert_called_once_with(self.single_cluster_node.alignment,
                                                         self.prg_builder.min_match_length)

    @patch.object(SingleClusterNode, SingleClusterNode._infer_if_this_node_should_have_no_child.__name__, return_value=True)
    def test___get_children___node_should_have_no_child(self, *uninteresting_mocks):
        self.subsetup()
        node = SingleClusterNode(1, self.alignment, None, self.prg_builder, False)

        expected = []
        actual = node._get_children()

        self.assertEqual(expected, actual)

    @patch.object(SingleClusterNode, SingleClusterNode._infer_if_this_node_should_have_no_child.__name__, return_value=False)
    @patch.object(SingleClusterNode, SingleClusterNode._infer_if_should_not_cluster.__name__, return_value=True)
    def test___get_children___single_SingleClusterNode_child(self, *uninteresting_mocks):
        self.subsetup()

        # we mock get_children() here so that we don't trigger the recursion tree building
        with patch.object(SingleClusterNode, SingleClusterNode._get_children.__name__):
            node = SingleClusterNode(1, self.alignment, None, self.prg_builder, False)
        node.all_intervals = [Interval(IntervalType.Match, 0, 0)]

        # now we mock SingleClusterNode so that we don't trigger the recursion tree building
        child_mock = Mock()
        with patch("make_prg.recursion_tree.SingleClusterNode", return_value=child_mock) as SingleClusterNode_mock:
            expected = [child_mock]
            actual = node._get_children()
            self.assertEqual(expected, actual)
            SingleClusterNode_mock.assert_called_once()

            # assert_called_once_with() does not work because can't check if MSAs are equal
            kwargs = SingleClusterNode_mock.call_args_list[0].kwargs
            self.assertTrue(first_dict_contained_in_second({
                    "nesting_level": node.nesting_level+1,
                    "parent": node,
                    "prg_builder": self.prg_builder,
                    "force_no_child": True
                }, kwargs))
            self.assertTrue(equal_msas(self.alignment[:,0:1], kwargs["alignment"]))

    @patch.object(SingleClusterNode, SingleClusterNode._infer_if_this_node_should_have_no_child.__name__, return_value=False)
    @patch.object(SingleClusterNode, SingleClusterNode._infer_if_should_not_cluster.__name__, return_value=False)
    def test___get_children___single_MultiClusterNode_child(self, *uninteresting_mocks):
        self.subsetup()

        # we mock get_children() here so that we don't trigger the recursion tree building
        with patch.object(SingleClusterNode, SingleClusterNode._get_children.__name__):
            node = SingleClusterNode(1, self.alignment, None, self.prg_builder, False)
        node.all_intervals = [Interval(IntervalType.Match, 0, 0)]

        # now we mock MultiClusterNode so that we don't trigger the recursion tree building
        child_mock = Mock()
        with patch("make_prg.recursion_tree.MultiClusterNode", return_value=child_mock) as MultiClusterNode_mock:
            expected = [child_mock]
            actual = node._get_children()
            self.assertEqual(expected, actual)
            MultiClusterNode_mock.assert_called_once()

            # assert_called_once_with() does not work because can't check if MSAs are equal
            kwargs = MultiClusterNode_mock.call_args_list[0].kwargs
            self.assertTrue(first_dict_contained_in_second({
                    "nesting_level": node.nesting_level+1,
                    "parent": node,
                    "prg_builder": self.prg_builder,
                    "force_no_child": False
                }, kwargs))
            self.assertTrue(equal_msas(self.alignment[:,0:1], kwargs["alignment"]))

    @patch.object(SingleClusterNode, SingleClusterNode._infer_if_this_node_should_have_no_child.__name__, return_value=False)
    @patch.object(SingleClusterNode, SingleClusterNode._infer_if_should_not_cluster.__name__, side_effect=[True, False, True, False])
    def test___get_children___multiple_children(self, *uninteresting_mocks):
        self.subsetup()

        # we mock get_children() here so that we don't trigger the recursion tree building
        with patch.object(SingleClusterNode, SingleClusterNode._get_children.__name__):
            node = SingleClusterNode(1, self.alignment, None, self.prg_builder, False)
        node.all_intervals = [Interval(IntervalType.Match, 0, 0),
                              Interval(IntervalType.NonMatch, 1, 1),
                              Interval(IntervalType.Match, 2, 2),
                              Interval(IntervalType.NonMatch, 3, 3)]

        # now we mock SingleClusterNode and MultiClusterNode so that we don't trigger the recursion tree building
        single_cluster_node_mock_1 = Mock()
        single_cluster_node_mock_2 = Mock()
        multi_cluster_node_mock_1 = Mock()
        multi_cluster_node_mock_2 = Mock()
        with patch("make_prg.recursion_tree.MultiClusterNode", side_effect=[multi_cluster_node_mock_1, multi_cluster_node_mock_2]), \
            patch("make_prg.recursion_tree.SingleClusterNode", side_effect=[single_cluster_node_mock_1, single_cluster_node_mock_2]):
            expected = [single_cluster_node_mock_1, multi_cluster_node_mock_1, single_cluster_node_mock_2, multi_cluster_node_mock_2]
            actual = node._get_children()
            self.assertEqual(expected, actual)
