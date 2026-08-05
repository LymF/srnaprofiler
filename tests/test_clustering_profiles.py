"""The merge needs a profile per island, computed from the reads in it."""
from features import island_profiles
from samio import Aln


def _reads(seq, positions, strand="+"):
    return [Aln("c", p, strand, seq, 1) for p in positions]


def test_profile_reports_dominant_size_and_strand():
    alns = _reads("A" * 24, [100, 110, 120])
    prof = island_profiles([(100, 150)], alns)
    assert prof[0]["size"] == 24
    assert prof[0]["frac_top"] == 1.0


def test_profile_of_dual_strand_island():
    alns = _reads("A" * 24, [100, 110]) + _reads("A" * 24, [120, 130], strand="-")
    prof = island_profiles([(100, 160)], alns)
    assert prof[0]["frac_top"] == 0.5


def test_island_without_reads_has_unknown_profile():
    prof = island_profiles([(5000, 5100)], _reads("A" * 24, [100]))
    assert prof[0]["size"] is None
    assert prof[0]["frac_top"] is None


def test_dominant_size_wins_over_minority():
    alns = _reads("A" * 24, [100, 110, 120]) + _reads("A" * 21, [130])
    prof = island_profiles([(100, 160)], alns)
    assert prof[0]["size"] == 24
