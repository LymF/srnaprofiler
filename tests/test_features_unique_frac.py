"""How much of a locus's read mass comes from uniquely-mapping reads.

57% of the piRNA loci have none, so their whole signature is borrowed from
repeats. multimap_frac counts alignment LINES, which over-weights the very
reads that are least certain; this one is on the same 1/NH scale as `reads`.
"""
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _r(nh, n=1, start=100):
    return [Aln("c", start + i, "+", "A" * 24, nh) for i in range(n)]


def test_all_unique_gives_one():
    r = build_record(("c", 0, 300), _r(1, 5), 1_000_000, CFG)
    assert r.unique_frac == 1.0


def test_all_multimapping_gives_zero():
    r = build_record(("c", 0, 300), _r(5, 5), 1_000_000, CFG)
    assert r.unique_frac == 0.0


def test_it_is_weighted_not_a_line_count():
    # 1 unique read + 4 lines of 4x multimappers: weighted mass is 1 + 4*(1/4) = 2,
    # so half of it is unique -- while the line count would say 1/5 = 0.2
    alns = _r(1, 1) + _r(4, 4, start=200)
    r = build_record(("c", 0, 300), alns, 1_000_000, CFG)
    assert r.unique_frac == 0.5
    assert r.multimap_frac == 0.8          # the raw line-based diagnostic is unchanged


def test_empty_locus_is_zero():
    assert build_record(("c", 0, 300), [], 1_000_000, CFG).unique_frac == 0.0
