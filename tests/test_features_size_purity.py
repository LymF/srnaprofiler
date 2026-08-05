"""A locus can hold two size populations; report it instead of hiding it.

18% of loci have dominant-size purity below 0.5, already at pad=75, so this is
biology we were discarding by reporting only the dominant DicerCall.
"""
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _reads(length, n, start=100):
    return [Aln("c", start + i, "+", "A" * length, 1) for i in range(n)]


def test_pure_locus_reports_purity_one_and_no_second_size():
    r = build_record(("c", 100, 300), _reads(24, 50), 1_000_000, CFG)
    assert r.size_purity == 1.0
    assert r.dicer_call_2 == ""


def test_bimodal_locus_reports_both():
    r = build_record(("c", 100, 300), _reads(24, 60) + _reads(21, 40, start=200),
                     1_000_000, CFG)
    assert 0.55 < r.size_purity < 0.65
    assert r.dicer_call == "24"
    assert r.dicer_call_2 == "21"


def test_second_size_omitted_when_purity_is_high():
    r = build_record(("c", 100, 300), _reads(24, 95) + _reads(21, 5, start=200),
                     1_000_000, CFG)
    assert r.size_purity > 0.9
    assert r.dicer_call_2 == ""


def test_empty_locus_is_safe():
    r = build_record(("c", 100, 300), [], 1_000_000, CFG)
    assert r.size_purity == 0.0
    assert r.dicer_call_2 == ""
