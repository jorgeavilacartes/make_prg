from unittest import TestCase
from unittest.mock import patch
from tests.test_helpers import make_alignment, equal_msas
from make_prg.recursion_tree import MultiClusterNode, SingleClusterNode
from make_prg.prg_builder import PrgBuilder
from make_prg.from_msa import MSA
from pathlib import Path

@patch.object(PrgBuilder, PrgBuilder.get_next_node_id.__name__, return_value=0)
@patch.object(SingleClusterNode, PrgBuilder.__init__.__name__, return_value=None)
@patch("make_prg.prg_builder.load_alignment_file")
class TestMultiClusterNode(TestCase):
    # can't apply patches to self.setUp()...
    def setup_get_sub_alignment_by_list_id(self) -> None:
        self.alignment = make_alignment(
            ["AAAT", "C--C", "AATT", "GNGG"], ["s1", "s2", "s3", "s4"]
        )
        self.prg_builder = PrgBuilder("locus", Path("msa"), "fasta", 5, 7)
        self.multi_cluster_node = MultiClusterNode(1, self.alignment, None, self.prg_builder, False)

    def test___get_sub_alignment_by_list_id___GivenOrderedIds_SubalignmentInSequenceOrder(self, *mocks):
        self.setup_get_sub_alignment_by_list_id()
        expected = MSA([self.alignment[0], self.alignment[2]])
        actual = self.multi_cluster_node._get_sub_alignment_by_list_id(["s1", "s3"])
        self.assertTrue(equal_msas(expected, actual))

    def test___get_sub_alignment_by_list_id___GivenUnorderedIds_SubalignmentStillInSequenceOrder(self, *mocks):
        """
        Sequences given rearranged are still output in input order
        """
        self.setup_get_sub_alignment_by_list_id()
        expected = MSA([self.alignment[0], self.alignment[2]])
        actual = self.multi_cluster_node._get_sub_alignment_by_list_id(["s3", "s1"])
        self.assertTrue(equal_msas(expected, actual))
