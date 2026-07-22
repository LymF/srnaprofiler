"""pingpong_frac: bounded [0,1] fraction of opposite-strand pairs overlapping by 10 nt.

Fixture: 4 plus reads (5' at 100) and 4 minus reads (5' at 109) form 4x4 = 16
weighted pairs, all with a 5'-5' overlap of exactly 10 -> fraction 1.0. The raw z
is guarded to 0 here (no background), which is why the classifier uses the fraction.
"""
from samio import read_alignments
from features import build_record

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def test_pingpong_frac_computed():
    alns = list(read_alignments("tests/data/pingpong_strong.sam"))
    r = build_record(("chr1", 50, 200), alns, lib_size=None, config=CFG)
    assert r.pingpong_pairs10 == 16     # 4 x 4 weighted pairs at overlap 10
    assert r.pingpong_frac == 1.0       # every pair overlaps by 10
    assert r.pingpong_z == 0.0          # guarded (no background) -> diagnostic only
