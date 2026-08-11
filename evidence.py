"""Null tests for the per-locus evidence scores."""
import math

from scipy import stats

# log10(1/p) of the smallest p we report; below this the test is saturated.
_MAX_SCORE = 300.0
# Fallback null for the 5' base when no genome is given. NOT a standard: proTRAC
# applies fraction thresholds (-1Tor10A 0.75) and runs no such test. With -g the
# null is the locus's own DNA composition instead (samio.base_null).
BASE_P0 = 0.25


def score_from_p(p):
    """log10(1/p), floored at 0 and capped at _MAX_SCORE."""
    if p is None or p >= 1.0:
        return 0.0
    if p <= 0.0:
        return _MAX_SCORE
    # -log10(p), not log10(1/p): the reciprocal overflows for subnormal p.
    return round(min(-math.log10(p), _MAX_SCORE), 3)


def effect_strength(k, n, p0, conf=0.95):
    """How far the effect clears its null, in [0, 1]: Wilson lower bound of k/n at
    `conf`, rescaled as (L - p0) / (1 - p0). Does not grow with n, so the same
    biology scores the same in a small locus and a large one.
    """
    if p0 is None or n is None:
        return 0.0
    n = int(round(n))
    k = int(round(k))
    if n <= 0 or k <= 0 or k > n or not (0.0 <= p0 < 1.0):
        return 0.0
    z = stats.norm.isf((1.0 - conf) / 2.0)
    phat = k / n
    denom = 1.0 + z * z / n
    centre = (phat + z * z / (2 * n)) / denom
    half = z * math.sqrt(phat * (1 - phat) / n + z * z / (4 * n * n)) / denom
    lower = centre - half
    return round(max(0.0, (lower - p0) / (1.0 - p0)), 4)


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


def strand_bidirectional_score(frac_top, n, null=0.10):
    """Evidence that BOTH strands produce: minority strand above `null`.

    One-sided binomial. Bidirectional does not mean 50/50 -- a cluster can be
    genuinely bidirectional at 70/30 -- so this asks whether the minority strand
    carries more than a uni-strand locus would leak, not whether the split is even.
    """
    n = int(round(n))
    if n <= 0 or frac_top is None:
        return 0.0
    k = int(round(min(frac_top, 1.0 - frac_top) * n))
    if k <= 0:
        return 0.0
    return score_from_p(stats.binomtest(k, n, null, alternative="greater").pvalue)


def strand_balance_score(frac_top, n, delta=0.15):
    """Equivalence (TOST) that the true strand split lies within 0.5 +/- delta.

    Two one-sided tests: H0a p>=0.5+delta ("less"), H0b p<=0.5-delta ("greater").
    TOST takes the LARGER p-value, so both one-sided nulls must be rejected.
    """
    n = int(round(n))
    if n <= 0 or frac_top is None:
        return 0.0
    k = int(round(frac_top * n))
    if k < 0 or k > n:
        return 0.0
    p_hi = stats.binomtest(k, n, 0.5 + delta, alternative="less").pvalue
    p_lo = stats.binomtest(k, n, 0.5 - delta, alternative="greater").pvalue
    return score_from_p(max(p_hi, p_lo))


def pingpong_score_empirical(counts, target=10, lo=1, hi=20):
    """Excess at the 10-nt offset over the rate of the locus's own other offsets,
    estimated from lo..hi excluding `target`. A uniform 1/20 null would read
    coverage structure as ping-pong.
    """
    if not counts:
        return 0.0
    bg = [counts.get(o, 0.0) for o in range(lo, hi + 1) if o != target]
    peak = counts.get(target, 0.0)
    if peak <= 0 or not bg:
        return 0.0
    rate = sum(bg) / len(bg)
    if rate <= 0:
        # No background at all: one offset carries everything. Real, but a rate
        # cannot be estimated, so report nothing rather than infinity.
        return 0.0
    return score_from_p(stats.poisson.sf(int(round(peak)) - 1, rate))


def pingpong_score(pairs_at_10, total_pairs, n_offsets=20):
    """Excess of 5'-5' pairs at the 10-nt offset over the other offsets."""
    k, n = int(round(pairs_at_10)), int(round(total_pairs))
    if n <= 0 or k <= 0 or k > n:
        return 0.0
    return score_from_p(
        stats.binomtest(k, n, 1.0 / n_offsets, alternative="greater").pvalue)


def density_score(reads, span_bp, lib_reads, genome_bp):
    # NOT wired into the pipeline: it restated the read count and was high by
    # construction, since a locus only exists after clearing --mincov.
    """Read pile-up over what a uniform distribution would put in this span."""
    if span_bp <= 0 or lib_reads <= 0 or genome_bp <= 0 or reads <= 0:
        return 0.0
    expected = lib_reads * (span_bp / genome_bp)
    if expected <= 0:
        return 0.0
    return score_from_p(stats.poisson.sf(int(round(reads)) - 1, expected))
