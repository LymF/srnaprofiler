from samio import read_alignments
from features import build_record, revcomp

CONFIG = {"dicermin": 20, "dicermax": 24, "strand_cutoff": 0.8}


def test_revcomp():
    assert revcomp("TAGC") == "GCTA"


def test_build_record_core():
    alns = list(read_alignments("tests/data/toy.sam"))
    region = ("chr1", 100, 160)
    r = build_record(region, alns, lib_size=1_000_000, config=CONFIG)
    assert r.reads == 5
    assert r.frac_top == 1.0
    assert r.strand_call == "+"
    assert r.dicer_call == "21"          # 4/5 = 0.8 in the window
    assert r.short == 0 and r.long == 1  # the 30nt one is "long"
    assert round(r.frac_in_window, 2) == 0.8
    assert r.pct_1u == round(4 / 5, 4)    # r1..r4 start with T(=U); r5 with G
    assert r.reads_rpm == 5.0             # 5/1e6*1e6
    assert r.multimap_frac == round(1 / 5, 4)  # r4 has NH=2
