"""Island depth is on the 1/NH scale, like every other count in the pipeline.

The threshold comes from the weighted library size, so comparing it against raw
alignment lines made it ~3x looser than declared on a bowtie1 -k run.
"""
from features import chromosome_depth
from samio import Aln


def test_unique_reads_count_one_each():
    d = chromosome_depth([Aln("c", 10, "+", "AAAAA", 1)])
    assert d[10] == 1.0
    assert d[14] == 1.0
    assert 15 not in d


def test_multimapper_counts_a_fraction():
    d = chromosome_depth([Aln("c", 10, "+", "AAAAA", 4)])
    assert d[10] == 0.25


def test_four_copies_of_a_4x_multimapper_make_one():
    d = chromosome_depth([Aln("c", 10, "+", "AAAAA", 4)] * 4)
    assert d[10] == 1.0


def test_depth_accumulates_across_reads():
    d = chromosome_depth([Aln("c", 10, "+", "AAAAA", 1),
                          Aln("c", 12, "+", "AAAAA", 2)])
    assert d[12] == 1.5
