# tests/test_features_signals.py
from samio import read_alignments
from features import build_record

CONFIG = {"dicermin": 8, "dicermax": 12, "strand_cutoff": 0.8}


def test_pingpong_pair_detected():
    alns = list(read_alignments("tests/data/pingpong.sam"))
    r = build_record(("chr1", 100, 160), alns, lib_size=None, config=CONFIG)
    assert r.pingpong_pairs10 == 1     # one +/- pair with overlap 10
