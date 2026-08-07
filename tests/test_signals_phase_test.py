"""Windowed phasing test (PhasiHunter convention: 11 cycles).

The p-value, the count of in-phase cycles and the population must all come from
the SAME window. Mixing a locus-wide count with a windowed population is what
produced p-values of exactly 0, which a hypergeometric can never yield.
"""
import pytest
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


def test_scan_does_not_prefer_an_empty_window():
    # A real phased run plus a lone read far away. Ranking windows by p-value
    # picked the window holding just that read -- P(X>=1) with one draw beats the
    # p of a window with a genuine seven-cycle run -- and reported number=1.
    positions = _phased(7, start=1000) + [(90000, 2.0)]
    res = phase_test(positions, region_start=1000, phase_length=21)
    assert res.number == 7
    assert res.score > 1.0


"""The scan tests one window per cycle and reports the best, so the p-value of
that window must pay for the selection. Uncorrected, the false-positive rate
grew with locus length: 2x inflation at ~1 window, 33x at ~134."""


def _in_register(start, n_cycles, phase=21, per_cycle=4):
    """Reads in perfect register, one stack per cycle."""
    return [(start + c * phase, float(per_cycle)) for c in range(n_cycles)]


def test_single_window_locus_is_not_corrected():
    # span fits in one window: nothing was selected, nothing to pay for
    pos = _in_register(1000, 8)
    r = phase_test(pos, region_start=1000, phase_length=21, cycles=11)
    assert r.pvalue == pytest.approx(_raw_pvalue(pos, 1000, 21, 8), rel=1e-6)


def test_scanned_locus_pays_for_the_windows_it_tried():
    # 60 cycles of span with the phased block at the start: the scan tries ~60
    # windows, so the reported p must be worse than the raw window p
    pos = _in_register(1000, 11) + [(1000 + 40 * 21 + 7, 1.0), (1000 + 59 * 21 + 3, 1.0)]
    r = phase_test(pos, region_start=1000, phase_length=21, cycles=11)
    raw = _raw_pvalue(_in_register(1000, 11), 1000, 21, 11)
    assert r.pvalue > raw
    assert r.pvalue <= 1.0


def test_corrected_p_is_capped_at_one():
    pos = [(1000 + i * 7, 1.0) for i in range(200)]     # noise over a long span
    r = phase_test(pos, region_start=1000, phase_length=21, cycles=11)
    assert r.pvalue <= 1.0


def _raw_pvalue(positions, start, phase, cycles):
    from signals import _phase_window
    return _phase_window(sorted(positions), start, phase, cycles).pvalue
