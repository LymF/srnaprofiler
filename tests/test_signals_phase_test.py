"""Windowed phasing test (PhasiHunter convention: 11 cycles).

The p-value, the count of in-phase cycles and the population must all come from
the SAME window. Mixing a locus-wide count with a windowed population is what
produced p-values of exactly 0, which a hypergeometric can never yield.
"""
from signals import phase_test


def _phased(n, start=1000, step=21, abundance=10.0):
    return [(start + i * step, abundance) for i in range(n)]


def test_perfectly_phased_locus_is_significant():
    res = phase_test(_phased(11), region_start=1000, phase_length=21)
    assert res.number == 11
    assert res.ratio == 1.0
    assert res.pvalue < 0.01


def test_scattered_positions_are_not_significant():
    # 11 reads spread over the same span but off-register
    positions = [(1000 + i * 21 + (i % 7), 10.0) for i in range(11)]
    res = phase_test(positions, region_start=1000, phase_length=21)
    assert res.pvalue > 0.01


def test_pvalue_is_never_exactly_zero():
    # a locus far longer than the window: the in-phase count must be clipped to
    # the window, so the hypergeometric stays well defined
    res = phase_test(_phased(60), region_start=1000, phase_length=21)
    assert res.pvalue > 0.0
    assert res.number <= 11


def test_window_shorter_than_11_cycles():
    res = phase_test(_phased(3), region_start=1000, phase_length=21)
    assert res.number == 3
    assert 0.0 < res.pvalue <= 1.0


def test_empty_locus():
    res = phase_test([], region_start=0, phase_length=21)
    assert (res.score, res.number, res.ratio, res.pvalue) == (0.0, 0, 0.0, 1.0)


def test_best_window_is_picked():
    # noise, then a clean phased run further along: the scan must find the run
    noise = [(1000 + i * 5, 1.0) for i in range(10)]
    res = phase_test(noise + _phased(11, start=5000), region_start=1000, phase_length=21)
    assert res.number == 11
    assert res.pvalue < 0.01
