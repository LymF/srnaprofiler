from samio import Aln
from class_signals import coverage_track


def test_coverage_track_strand_lanes_and_clip():
    alns = [Aln("chr", 0, "+", "A" * 50, 1), Aln("chr", 60, "-", "A" * 20, 1)]
    sense, anti = coverage_track(alns, 0, 100)
    assert len(sense) == 100 and len(anti) == 100
    assert sense[0] == 1 and sense[49] == 1 and sense[50] == 0
    assert anti[60] == 1 and anti[79] == 1 and anti[80] == 0
    assert all(v == 0 for v in anti[:60])


def test_coverage_track_is_raw_and_stacks():
    # nh=2 still counts +1 (raw, not 1/NH); two overlapping + reads stack to 2
    alns = [Aln("chr", 0, "+", "A" * 10, 2), Aln("chr", 0, "+", "A" * 10, 1)]
    sense, _ = coverage_track(alns, 0, 10)
    assert sense[0] == 2


def test_coverage_track_empty_region():
    assert coverage_track([], 5, 5) == ([], [])
