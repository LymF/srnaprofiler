"""Null tests for the per-locus evidence scores."""
import math

from scipy import stats

# log10(1/p) of the smallest p we report; below this the test is saturated.
_MAX_SCORE = 300.0
# proTRAC tests the 5' base against random base composition.
BASE_P0 = 0.25


def score_from_p(p):
    """log10(1/p), floored at 0 and capped at _MAX_SCORE."""
    if p is None or p >= 1.0:
        return 0.0
    if p <= 0.0:
        return _MAX_SCORE
    return round(min(math.log10(1.0 / p), _MAX_SCORE), 3)


def base_bias_score(k, n, p0=BASE_P0):
    """Enrichment of a 5' base over chance (one-sided binomial)."""
    k, n = int(round(k)), int(round(n))
    if n <= 0 or k <= 0 or k > n:
        return 0.0
    return score_from_p(stats.binomtest(k, n, p0, alternative="greater").pvalue)


def strand_bias_score(frac_top, n):
    """Departure from an even strand split (two-sided binomial)."""
    n = int(round(n))
    if n <= 0 or frac_top is None:
        return 0.0
    k = int(round(frac_top * n))
    if k < 0 or k > n:
        return 0.0
    return score_from_p(stats.binomtest(k, n, 0.5, alternative="two-sided").pvalue)
