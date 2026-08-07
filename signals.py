"""Biogenesis signals: ping-pong z-score and phasing."""
import bisect
import math
import statistics
from collections import Counter, defaultdict, namedtuple
from math import comb


def offset_spectrum(a5, b5, lo=1, hi=20):
    """5'-5' offset spectrum: Counter[d] = sum over p of a5[p] * b5[p + d - 1].

    `d` is the overlap length in `lo..hi` (ping-pong peak at d=10). Cross-strand
    ping-pong passes plus5/minus5; same-strand phasing passes the same dict twice,
    which yields the 5'-5' distance distribution (peaks at multiples of the phase
    length) and is phase-invariant, so it pools correctly across loci.
    """
    counts = Counter()
    for p, ac in a5.items():
        for d in range(lo, hi + 1):
            bc = b5.get(p + d - 1, 0)
            if bc:
                counts[d] += ac * bc
    return counts


def pingpong_zscore(counts, target=10, lo=1, hi=20, min_bg=3, cap=100.0):
    """Z-score of the 5'-5' overlap at `target` (default 10 nt) against the background.

    counts: {overlap_length: n_pairs}. Background = offsets in [lo,hi] minus the target.

    Guarded and capped: on a narrow, ultra-dense locus the background collapses to a
    couple of non-zero offsets, its stdev goes to ~0 and the raw z explodes into the
    thousands, which is not comparable across loci. We require at least `min_bg`
    non-zero background offsets for the z to be meaningful, and cap its magnitude at
    `cap`. (For classification use the bounded `pingpong_frac` instead — this z is a
    diagnostic.)
    """
    background = [counts.get(o, 0) for o in range(lo, hi + 1) if o != target]
    peak = counts.get(target, 0)
    nonzero_bg = sum(1 for v in background if v > 0)
    if len(background) < 2 or nonzero_bg < min_bg:
        return 0.0, peak
    sd = statistics.pstdev(background)
    if sd == 0:
        return 0.0, peak
    mean = statistics.fmean(background)
    z = (peak - mean) / sd
    z = max(-cap, min(cap, z))
    return round(z, 3), peak


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


PhaseResult = namedtuple("PhaseResult", ["score", "number", "ratio", "pvalue"])


def phase_test(positions, region_start, phase_length, cycles=11):
    """Phasing over the best window of `cycles` phase cycles.

    positions: [(pos5, abundance)], 5' ends already pooled across strands.

    Scans windows of `cycles * phase_length` nt in steps of one cycle and keeps
    the most significant one, so local phasing is not diluted by a long island.
    Everything the hypergeometric consumes comes from that one window:

      population = 2 * window_len - 1   positions a read could start at
      successes  = 2 * cycles - 1       positions that are in register
      draws      = distinct positions occupied in the window
      observed   = of those, how many are in the dominant register

    The 2n-1 convention is PhasiHunter's (both strands can open a read at a
    position); what matters is that successes/population ~ 1/phase_length, the
    chance a random position lands in register.
    """
    if not positions:
        return PhaseResult(0.0, 0, 0.0, 1.0)

    span = max(p for p, _ in positions) - min(p for p, _ in positions) + 1
    window_len = cycles * phase_length
    if span <= window_len:
        return _phase_window(positions, region_start, phase_length,
                             min(cycles, max(1, -(-span // phase_length))))
    # Sort once and slice each window by binary search: refiltering the whole list
    # per window is O(windows x positions), which is untenable on a large island.
    ordered = sorted(positions)
    coords = [p for p, _ in ordered]
    lo, hi = coords[0], coords[-1]
    best = None
    tried = 0
    for start in range(lo, hi + 1, phase_length):
        i = bisect.bisect_left(coords, start)
        j = bisect.bisect_left(coords, start + window_len)
        if i == j:
            continue
        tried += 1
        res = _phase_window(ordered[i:j], region_start, phase_length, cycles)
        # Rank by evidence, never by p-value: the hypergeometric rewards
        # near-empty windows, so ranking by p picked the emptiest one.
        if best is None or (res.number, res.score) > (best.number, best.score):
            best = res
    if best is None:
        return PhaseResult(0.0, 0, 0.0, 1.0)
    # The reported window was the best of `tried`, so its p must pay for the
    # selection. Bonferroni is conservative here (the windows overlap and are
    # correlated), which is the right direction for evidence.
    return best._replace(pvalue=min(1.0, best.pvalue * tried))


def _phase_window(positions, region_start, phase_length, cycles):
    """Phasing metrics and p-value for a single window."""
    score, number, ratio = phase_metrics(positions, region_start, phase_length)
    window_len = cycles * phase_length
    pvalue = hypergeometric_phase_pvalue(
        population=max(2 * window_len - 1, 1),
        phased_positions=max(2 * cycles - 1, 1),
        draws=len({p for p, _ in positions}),
        observed=min(number, max(2 * cycles - 1, 1)))
    return PhaseResult(score, number, ratio, pvalue)


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
