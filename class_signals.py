"""Per-class signal accumulators over samio.Aln (1/NH-weighted, RNA 5'->3')."""
import statistics
from collections import Counter
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


def pattern_row(alns, minlen=15, maxlen=35):
    """42 z-scored size x strand counts (0-20 = +, 21-41 = -), weighted by 1/NH."""
    span = maxlen - minlen + 1
    values = [0.0] * (span * 2)
    for a in alns:
        L = len(_oriented(a))
        if L < minlen or L > maxlen:
            continue
        idx = L - minlen if a.strand == "+" else (L - minlen) + span
        values[idx] += 1.0 / a.nh
    nonzero = [v for v in values if v != 0.0]
    if len(nonzero) < 2:
        return [0.0] * len(values)
    mean = statistics.fmean(values)
    sd = statistics.pstdev(values)
    if sd == 0:
        return [0.0] * len(values)
    return [(v - mean) / sd for v in values]


def coverage_profile(alns, region_start, region_end, bins=100):
    """Per-strand weighted coverage over the locus span, rescaled to `bins` bins.

    Each bin holds the mean per-position 1/NH weight within its slice of the
    region (sum of weights divided by positions-per-bin), so the result is
    comparable across loci of different lengths regardless of `bins`.
    """
    length = region_end - region_start
    prof = {"+": [0.0] * bins, "-": [0.0] * bins}
    if length <= 0:
        return prof
    for a in alns:
        w = 1.0 / a.nh
        lane = prof[a.strand]
        read_start = max(a.pos, region_start)
        read_end = min(a.pos + len(a.seq), region_end)
        for p in range(read_start, read_end):
            b = int((p - region_start) / length * bins)
            if b >= bins:
                b = bins - 1
            lane[b] += w
    positions_per_bin = length / bins
    for lane in prof.values():
        for b in range(bins):
            lane[b] /= positions_per_bin
    return prof


def strand_sequences(alns, strand, size_lo, size_hi, trim):
    """Distinct RNA-oriented, size-filtered, trimmed, T->U sequences on `strand`."""
    seen = set()
    for a in alns:
        if a.strand != strand:
            continue
        seq = _oriented(a)
        L = len(seq)
        if L < size_lo or L > size_hi:
            continue
        trimmed = seq[:trim].replace("T", "U").replace("t", "u")
        if len(trimmed) == trim:
            seen.add(trimmed)
    return sorted(seen)


def write_fasta(sequences, path):
    with open(path, "w") as fh:
        for i, seq in enumerate(sequences, 1):
            fh.write(f">seq_{i}\n{seq}\n")
