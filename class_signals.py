"""Per-class signal accumulators over samio.Aln (1/NH-weighted, RNA 5'->3')."""
from collections import Counter, defaultdict
from features import revcomp

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
