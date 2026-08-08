"""Size purity is a band, not a single length.

A Dicer product has a +/-1 nt spread, so scoring purity on one length alone
called a clean 21-nt locus impure: the median size_purity was 0.196 in
Drosophila and dicer_call_2 fired on 98.7% of the loci, half of them naming a
neighbouring length -- the same population's spread, not a second pathway.
"""
from records import ClusterRecord
from features import summarize_sizes


def _sizes(counts):
    return {str(k): float(v) for k, v in counts.items()}


def test_the_dominant_band_includes_the_neighbours():
    r = ClusterRecord(chrom="c", start=0, end=100)
    summarize_sizes(r, _sizes({20: 10, 21: 70, 22: 20}))
    assert r.size_purity == 1.0


def test_a_second_population_is_not_part_of_the_band():
    r = ClusterRecord(chrom="c", start=0, end=100)
    summarize_sizes(r, _sizes({21: 60, 24: 40}))
    assert r.size_purity == 0.6
    assert r.dicer_call_2 == "24"


def test_a_neighbouring_length_is_not_a_second_call():
    r = ClusterRecord(chrom="c", start=0, end=100)
    summarize_sizes(r, _sizes({21: 60, 22: 40}))
    assert r.dicer_call_2 == ""


def test_a_minor_second_peak_is_not_reported():
    """Below a fifth of the mass it is background, not a population."""
    r = ClusterRecord(chrom="c", start=0, end=100)
    summarize_sizes(r, _sizes({21: 90, 28: 10}))
    assert r.dicer_call_2 == ""


def test_the_second_call_is_the_peak_of_the_other_band():
    """Not merely the runner-up length: the top of the second population."""
    r = ClusterRecord(chrom="c", start=0, end=100)
    summarize_sizes(r, _sizes({21: 40, 22: 25, 27: 20, 28: 30}))
    assert r.dicer_call_2 == "28"


def test_a_single_length_locus_is_pure():
    r = ClusterRecord(chrom="c", start=0, end=100)
    summarize_sizes(r, _sizes({24: 50}))
    assert r.size_purity == 1.0
    assert r.dicer_call_2 == ""


def test_no_reads_in_the_window_leaves_the_columns_empty():
    r = ClusterRecord(chrom="c", start=0, end=100)
    summarize_sizes(r, {})
    assert r.size_purity == 0.0
    assert r.dicer_call_2 == ""
