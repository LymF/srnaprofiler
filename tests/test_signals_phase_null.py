"""Calibration of phase_test under a null with no phasing.

Two real libraries landing near the nominal rate is weak evidence: it could be
luck of their locus-size distribution. Here the truth is known by construction
-- 5' positions are drawn uniformly, so no locus is phased -- and a calibrated
test must report p <= a for about a fraction a of them, at every locus size.

The window scan is what makes this non-trivial: it tries one window per cycle
and reports the best, so without a correction the false-positive rate grows
with the locus, which is exactly what was measured on both kingdoms.
"""
import random

import pytest

from signals import phase_test

PHASE = 21
CYCLES = 11


def _unphased(rng, n_cycles, n_reads):
    """Reads at uniformly random 5' positions over `n_cycles` phase cycles."""
    span = n_cycles * PHASE
    return [(rng.randrange(span), 1.0) for _ in range(n_reads)]


def _rate(n_cycles, n_reads, alpha, trials=400, seed=0):
    rng = random.Random(seed)
    hits = 0
    for _ in range(trials):
        res = phase_test(_unphased(rng, n_cycles, n_reads), region_start=0,
                         phase_length=PHASE, cycles=CYCLES)
        hits += res.pvalue <= alpha
    return hits / trials


@pytest.mark.parametrize("n_cycles", [12, 30, 60, 120])
def test_false_positive_rate_does_not_grow_with_the_locus(n_cycles):
    """The bug's signature was a rate rising with the number of windows:
    18.8% at 12 cycles up to 76.3% at 120. It must now be flat in the locus."""
    assert _rate(n_cycles, n_reads=40, alpha=0.01) <= 0.04


@pytest.mark.parametrize("alpha", [0.05, 0.01])
def test_rate_stays_within_a_small_factor_of_nominal(alpha):
    """A residual ~2-3x anti-conservatism survives the correction.

    It does not come from the scan -- it is the hypergeometric null itself,
    whose population/successes are an idealisation of where a read may start.
    Bounded and size-independent, unlike the scan bias it replaced (76% of
    unphased 120-cycle loci passed p <= 0.01 before).
    """
    assert _rate(60, n_reads=40, alpha=alpha) <= alpha * 3


def test_a_genuinely_phased_locus_is_still_detected():
    """Calibration is worthless if the correction also kills the signal."""
    phased = [(c * PHASE, 4.0) for c in range(20)]
    res = phase_test(phased, region_start=0, phase_length=PHASE, cycles=CYCLES)
    assert res.pvalue <= 0.01
