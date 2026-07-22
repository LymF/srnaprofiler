"""Biogenesis signals: ping-pong z-score and phasing."""
import math
import statistics
from collections import defaultdict
from math import comb


def pingpong_zscore(counts, target=10, lo=1, hi=20):
    """Z-score of the 5'-5' overlap at `target` (default 10 nt) against the background.

    counts: {overlap_length: n_pairs}. Background = offsets in [lo,hi] minus the target.
    """
    background = [counts.get(o, 0) for o in range(lo, hi + 1) if o != target]
    peak = counts.get(target, 0)
    if len(background) < 2:
        return 0.0, peak
    sd = statistics.pstdev(background)
    if sd == 0:
        return 0.0, peak
    mean = statistics.fmean(background)
    return round((peak - mean) / sd, 3), peak


def phase_metrics(positions, region_start, phase_length):
    """PhasiHunter-style phase score: ratio * number * log(abun) at the dominant register.

    positions: [(pos5, abundance)]. Register = (pos5 - region_start) % phase_length.
    """
    reg_abun = defaultdict(int)
    reg_cycles = defaultdict(set)
    total = 0
    for pos5, ab in positions:
        r = (pos5 - region_start) % phase_length
        cyc = (pos5 - region_start) // phase_length
        reg_abun[r] += ab
        reg_cycles[r].add(cyc)
        total += ab
    if total == 0:
        return 0.0, 0, 0.0
    dom = max(reg_abun, key=reg_abun.get)
    phase_abun = reg_abun[dom]
    phase_number = len(reg_cycles[dom])
    phase_ratio = phase_abun / total
    score = phase_ratio * phase_number * math.log(phase_abun) if phase_abun > 0 else 0.0
    return round(score, 3), phase_number, round(phase_ratio, 3)


def hypergeometric_phase_pvalue(population, phased_positions, draws, observed):
    """P(X >= observed) drawing `draws` from `population` with `phased_positions` successes.

    comb() raises ValueError for negative arguments; the `min(phased_positions,
    draws)` upper bound on the loop ensures we only sum valid terms
    (i never exceeds the number of successes nor the number of draws).
    """
    # Degenerate cases (tiny clusters) make the hypergeometric ill-defined; a
    # cluster that cannot even fit the phased positions is never significant.
    if (observed <= 0 or population <= 0 or phased_positions <= 0
            or draws <= 0 or phased_positions > population or draws > population):
        return 1.0
    denom = comb(population, draws)
    if denom == 0:
        return 1.0
    p = 0.0
    for i in range(observed, min(phased_positions, draws) + 1):
        p += comb(phased_positions, i) * comb(population - phased_positions, draws - i) / denom
    return min(p, 1.0)
