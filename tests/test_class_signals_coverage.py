from samio import Aln
from class_signals import coverage_profile


def test_coverage_profile_bins_and_strand():
    # a + read covering the first half of a length-100 region
    alns = [Aln("chr", 0, "+", "A" * 50, 1)]
    prof = coverage_profile(alns, region_start=0, region_end=100, bins=10)
    assert len(prof["+"]) == 10
    assert prof["+"][0] > 0.0        # first bin covered
    assert prof["+"][9] == 0.0       # last bin not covered
    assert all(v == 0.0 for v in prof["-"])


def test_coverage_profile_weighting():
    alns = [Aln("chr", 0, "+", "A" * 100, 2)]   # w = 0.5 over whole region
    prof = coverage_profile(alns, region_start=0, region_end=100, bins=10)
    assert all(abs(v - 0.5) < 1e-9 for v in prof["+"])
