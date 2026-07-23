"""Per-class signal accumulators over samio.Aln (1/NH-weighted, RNA 5'->3')."""
from collections import Counter, defaultdict
from features import revcomp
from signals import offset_spectrum

_BASES = ("A", "C", "G", "T")


def _oriented(a):
    """Sequence in the RNA 5'->3' orientation (revcomp on the - strand)."""
    return a.seq if a.strand == "+" else revcomp(a.seq)


def base_composition(alns, minlen=15, maxlen=35):
    """First-base composition by size and strand, weighted by 1/NH."""
    comp = {s: {"+": {b: 0.0 for b in _BASES}, "-": {b: 0.0 for b in _BASES}}
            for s in range(minlen, maxlen + 1)}
    for a in alns:
        seq = _oriented(a)
        L = len(seq)
        if L < minlen or L > maxlen or not seq:
            continue
        first = seq[0]
        if first in comp[L][a.strand]:
            comp[L][a.strand][first] += 1.0 / a.nh
    return comp


def five_prime_ends(alns):
    """Weighted 5' end positions: (plus5, minus5) Counters."""
    plus5, minus5 = Counter(), Counter()
    for a in alns:
        w = 1.0 / a.nh
        if a.strand == "+":
            plus5[a.pos] += w
        else:
            minus5[a.pos + len(a.seq) - 1] += w
    return plus5, minus5


def pingpong_offsets(alns, hi=30):
    """Cross-strand 5'-5' overlap spectrum (ping-pong peak at 10)."""
    plus5, minus5 = five_prime_ends(alns)
    return offset_spectrum(plus5, minus5, lo=1, hi=hi)


def phasing_series(alns, strand, maxdist=100):
    """Same-strand 5'-5' distance spectrum for `strand` (phase-invariant).

    offset_spectrum keys are 1-indexed overlap length (= distance + 1); we remap
    by -1 so the result is keyed by true 5'-5' distance, peaking at multiples of
    the phase length (aligning with the R reference lines).
    """
    plus5, minus5 = five_prime_ends(alns)
    ends = plus5 if strand == "+" else minus5
    raw = offset_spectrum(ends, ends, lo=2, hi=maxdist + 1)
    return Counter({k - 1: v for k, v in raw.items()})
