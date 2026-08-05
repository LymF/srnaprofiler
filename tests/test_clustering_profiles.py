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


# A profile estimated from a handful of reads is noise, not evidence: the raw
# island has ~5 weighted reads, so its modal size flips with one read either way.
def test_profile_below_the_read_minimum_is_unknown():
    prof = island_profiles([(100, 160)], _reads("A" * 24, [100, 110]), min_reads=5)
    assert prof[0]["size"] is None
    assert prof[0]["frac_top"] is None


def test_profile_at_the_minimum_is_kept():
    prof = island_profiles([(100, 160)], _reads("A" * 24, [100, 110, 120, 130, 140]),
                           min_reads=5)
    assert prof[0]["size"] == 24


def test_minimum_counts_weighted_reads():
    # 8 alignments of a read mapping to 4 places = 2 weighted reads
    alns = [Aln("c", 100 + i, "+", "A" * 24, 4) for i in range(8)]
    assert island_profiles([(100, 160)], alns, min_reads=5)[0]["size"] is None


def test_default_keeps_every_profile():
    prof = island_profiles([(100, 160)], _reads("A" * 24, [100]))
    assert prof[0]["size"] == 24


# Which class window holds most of the read mass is far more stable with few
# reads than the exact modal length, and it is what the merge actually needs.
_WIN = [(20, 22), (24, 31)]


def test_window_survives_a_noisy_mode():
    # 3 reads of 25, 2 of 26, 2 of 27: no majority length, but all piRNA-window
    alns = (_reads("A" * 25, [100, 110, 120]) + _reads("A" * 26, [130, 140])
            + _reads("A" * 27, [150, 160]))
    prof = island_profiles([(100, 200)], alns, size_windows=_WIN)
    assert prof[0]["window"] == 1


def test_window_separates_the_pathways():
    p21 = island_profiles([(100, 200)], _reads("A" * 21, [100, 110, 120]),
                          size_windows=_WIN)
    p25 = island_profiles([(100, 200)], _reads("A" * 25, [100, 110, 120]),
                          size_windows=_WIN)
    assert p21[0]["window"] == 0 and p25[0]["window"] == 1


def test_no_dominant_window_is_unknown():
    # half 21-nt, half 25-nt: genuinely mixed, so no call
    alns = _reads("A" * 21, [100, 110]) + _reads("A" * 25, [120, 130])
    assert island_profiles([(100, 200)], alns, size_windows=_WIN)[0]["window"] is None


def test_mass_outside_every_window_is_unknown():
    # 19 nt degradation product: belongs to no pathway, must not be a barrier
    assert island_profiles([(100, 200)], _reads("A" * 19, [100, 110, 120]),
                           size_windows=_WIN)[0]["window"] is None


def test_window_is_absent_without_windows():
    assert island_profiles([(100, 200)], _reads("A" * 25, [100]))[0]["window"] is None
