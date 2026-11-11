"""Per-cluster figure data (sizes, nt composition, coverage, ping-pong) from a BAM."""
from collections import Counter

import pysam

import class_signals
import samio
import utils

_NT_POS = 26
_PP_HI = 30


def _r(x):
    return round(x, 3)


def _bin_lane(lane, bin_size):
    return [float(sum(lane[i:i + bin_size])) for i in range(0, len(lane), bin_size)]


def cluster_signals(alns, start, end, size_lo=18, size_hi=35, max_bins=120):
    """Figure data for one cluster window [start, end); see the module docstring."""
    alns = list(alns)
    n = size_hi - size_lo + 1
    plus, minus = [0.0] * n, [0.0] * n
    by_size = Counter()
    for a in alns:
        L = len(a.seq)
        if size_lo <= L <= size_hi:
            w = 1.0 / a.nh
            (plus if a.strand == "+" else minus)[L - size_lo] += w
        by_size[L] += 1.0 / a.nh

    # nucleotide composition of the dominant-size reads
    nt = [{"A": 0, "C": 0, "G": 0, "U": 0} for _ in range(_NT_POS)]
    if by_size:
        dom = max(by_size, key=lambda s: (by_size[s], -s))
        tot = [0.0] * _NT_POS
        cnt = [{"A": 0.0, "C": 0.0, "G": 0.0, "U": 0.0} for _ in range(_NT_POS)]
        for a in alns:
            if len(a.seq) != dom:
                continue
            w = 1.0 / a.nh
            for i, b in enumerate(class_signals._oriented(a)[:_NT_POS]):
                b = "U" if b in "Tt" else b.upper()
                if b in cnt[i]:
                    cnt[i][b] += w
                    tot[i] += w
        for i in range(_NT_POS):
            if tot[i] > 0:
                nt[i] = {b: _r(v / tot[i]) for b, v in cnt[i].items()}

    sense, anti = class_signals.coverage_track(alns, start, end)
    bin_size = max(1, -(-len(sense) // max_bins))
    pp = class_signals.pingpong_offsets(alns, hi=_PP_HI)
    return {
        "sizes": {"lo": size_lo, "plus": [_r(v) for v in plus],
                  "minus": [_r(v) for v in minus]},
        "nt": nt,
        "coverage": {"start": start, "bin": bin_size,
                     "plus": _bin_lane(sense, bin_size),
                     "minus": _bin_lane(anti, bin_size)},
        "pingpong": [_r(pp.get(o, 0.0)) for o in range(1, _PP_HI + 1)],
        "n_alns": len(alns),
    }


def signals_for(bam_path, clusters, xm_as_nh=None):
    """One figure-data entry per cluster name; a failing region is zero-filled."""
    if xm_as_nh is None:
        xm_as_nh = samio.prefers_xm_as_nh(bam_path)
    out = {}
    failed, first = 0, None
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        refs = set(fh.references)
        for c in clusters:
            alns = []
            # the TSV parser yields floats; coordinates must be ints here
            c_start, c_end = int(c.start), int(c.end)
            if c.chrom in refs:
                try:
                    alns = list(samio.read_region_handle(
                        fh, c.chrom, c_start, c_end, xm_as_nh=xm_as_nh))
                except Exception as e:
                    failed += 1
                    first = first or f"{type(e).__name__}: {e}"
            out[c.name] = cluster_signals(alns, c_start, c_end)
    if failed:
        utils.log(f"report figures: {failed} cluster(s) failed to read from the BAM "
                  f"(zero-filled); first error: {first}", "WARN")
    return out
