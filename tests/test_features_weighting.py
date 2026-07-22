"""Multi-mapping reads are counted fractionally (1/NH).

A read mapping to N places contributes 1/N at each locus. The fixture has 3
unique reads (NH=1) plus one read mapping to 4 places (NH=4), so the locus'
effective read count is 3 + 1/4 = 3.25, while multimap_frac stays a raw
diagnostic (1 of 4 alignment lines is a multimapper).
"""
from samio import read_alignments
from features import build_record

CONFIG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def test_reads_weighted_by_inverse_nh():
    alns = list(read_alignments("tests/data/multimap.sam"))
    r = build_record(("chr1", 100, 160), alns, lib_size=1_000_000, config=CONFIG)
    assert r.reads == 3.25                 # 3 * 1 + 1 * (1/4)
    assert r.reads_rpm == round(3.25, 4)   # 3.25 / 1e6 * 1e6
    assert r.multimap_frac == 0.25         # 1 of 4 alignment lines has NH>1
    # the size-21 window count is also weighted (all four reads are 21 nt)
    assert r.size_counts[21] == 3.25
