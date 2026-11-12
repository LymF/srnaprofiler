"""Per-cluster figure data for the HTML report, read in one pass from the BAM.

Each cluster region is fetched once and aggregated into two Counters, one keyed
by (strand, pos, length, nh, 5' base) and one by distinct
sequence (all lengths, for the sequences table); no per-read
list is kept. Memory is bounded by the number of distinct keys plus distinct
sequences: at most the number of records, and small for redundant small-RNA
data. The sequence Counter is not capped (an evicted sequence could return and
change the top rows); coverage and 5' ends are binned, so they cost O(bins),
not O(cluster length). Numbers match `class_signals`.
"""
from collections import Counter
import pysam

import samio
import utils
from features import revcomp
from signals import offset_spectrum

_NT_POS = 26
_PP_HI = 30
_PHASE_MAX = 100
_SKIP = 0x4 | 0x100 | 0x800  # unmapped, secondary, supplementary (samio._keep)


def _r(x):
    return round(float(x), 3)


def _accumulate(fh, chrom, start, end, *, xm_as_nh):
    """One pass over the region: (records kept, key Counter, sequence Counter).

    keys: (strand, pos, length, nh, 5' base) -> record count.
    seqs: (is_reverse, stored sequence) -> 1/NH-weighted reads.
    """
    keys, seqs = Counter(), Counter()
    n = 0
    for rec in fh.fetch(chrom, start, end):
        if rec.flag & _SKIP:
            continue
        seq = rec.query_sequence
        if not seq:
            continue
        nh = samio._nh(rec, xm_as_nh)
        rev = rec.is_reverse
        if rev:
            keys["-", rec.reference_start, len(seq), nh, revcomp(seq[-1])] += 1
        else:
            keys["+", rec.reference_start, len(seq), nh, seq[0]] += 1
        seqs[rev, seq] += 1.0 / nh
        n += 1
    return n, keys, seqs


def _add_cov(sums, full, a, b, c, B):
    """Add depth c over window offsets [a, b) to per-bin base sums.

    Partial first/last bins get c * overlap directly; fully covered bins go
    through the difference array `full` (c * B each).
    """
    ba, bb = a // B, (b - 1) // B
    if ba == bb:
        sums[ba] += c * (b - a)
        return
    sums[ba] += c * ((ba + 1) * B - a)
    sums[bb] += c * (b - bb * B)
    if bb > ba + 1:
        full[ba + 1] += c * B
        full[bb] -= c * B


def _cov_means(sums, full, B, length):
    out, run = [], 0
    for i, v in enumerate(sums):
        run += full[i]
        out.append(_r((v + run) / min(B, length - i * B)))
    return out


def _build(n, keys, seqs, start, end, *, dominant=None, size_lo=18, size_hi=35,
           max_bins=400, top_seqs=10):
    """Figure data from the counters of `_accumulate`."""
    span = size_hi - size_lo + 1
    sizes = {b: {"+": [0.0] * span, "-": [0.0] * span} for b in "ACGU"}
    by_size = Counter()
    ends = {"+": Counter(), "-": Counter()}
    length = max(end - start, 0)
    B = max(1, -(-length // max_bins))  # bin size
    nbins = -(-length // B)
    sums = {s: [0] * nbins for s in "+-"}
    full = {s: [0] * (nbins + 1) for s in "+-"}
    for (s, pos, L, nh, b), c in keys.items():
        w = c / nh
        b = "U" if b == "T" else b  # as base_composition: uppercase ACGT only
        if size_lo <= L <= size_hi:
            by_size[L] += w
            if b in sizes:
                sizes[b][s][L - size_lo] += w
        ends[s][pos if s == "+" else pos + L - 1] += w
        lo, hi = max(pos, start), min(pos + L, end)  # clipped like coverage_track
        if lo < hi:
            _add_cov(sums[s], full[s], lo - start, hi - start, c, B)

    if dominant is None or not size_lo <= dominant <= size_hi:
        dominant = max(by_size, key=lambda k: (by_size[k], -k)) if by_size else None

    cov, end_bins = {"start": start, "bin": B}, {"bin": B}
    for s in "+-":
        cov[s] = _cov_means(sums[s], full[s], B, length)
        lane = [0.0] * nbins
        for p, w in ends[s].items():
            if start <= p < end:
                lane[(p - start) // B] += w
        end_bins[s] = [_r(v) for v in lane]

    pp = offset_spectrum(ends["+"], ends["-"], lo=1, hi=_PP_HI)
    phasing = {}
    for s in "+-":
        # offset_spectrum keys are distance + 1 (see class_signals.phasing_series)
        raw = offset_spectrum(ends[s], ends[s], lo=2, hi=_PHASE_MAX + 1)
        phasing[s] = [_r(raw.get(k + 1, 0.0)) for k in range(1, _PHASE_MAX + 1)]

    # distinct sequences, RNA-oriented, T as U
    rna = Counter()
    for (rev, seq), w in seqs.items():
        o = (revcomp(seq) if rev else seq).upper().replace("T", "U")
        rna[o, "-" if rev else "+"] += w
    # ties compared at 1e-9 so float summation order cannot reorder them
    rows = sorted(rna.items(), key=lambda kv: (-round(kv[1], 9), kv[0][0], kv[0][1]))
    top = [[o, len(o), s, _r(w), o[0]] for (o, s), w in rows[:top_seqs]]

    logo = {}
    for s in "+-":
        cnt = [dict.fromkeys("ACGU", 0.0) for _ in range(_NT_POS)]
        for (o, st), w in rna.items():
            if st != s or len(o) != dominant:
                continue
            for i, b in enumerate(o[:_NT_POS]):
                if b in "ACGU":
                    cnt[i][b] += w
        tots = [sum(p.values()) for p in cnt]
        logo[s] = [{b: _r(v / t) for b, v in p.items()} if t > 0 else dict.fromkeys("ACGU", 0)
                   for p, t in zip(cnt, tots)]

    return {
        "n_alns": n,
        "dominant": dominant,
        "sizes": dict({"lo": size_lo},
                      **{b: {s: [_r(v) for v in sizes[b][s]] for s in "+-"} for b in "ACGU"}),
        "phasing": phasing,
        "pingpong": [_r(pp.get(o, 0.0)) for o in range(1, _PP_HI + 1)],
        "cov": cov,
        "ends": end_bins,
        "logo": logo,
        "seqs": top,
    }


def cluster_signals(fh, chrom, start, end, *, xm_as_nh, dominant=None, size_lo=18,
                    size_hi=35, max_bins=400, top_seqs=10):
    """Figure data for one cluster window [start, end) of an open AlignmentFile.

    `dominant` is the cluster's size call; without it (or outside the size
    window) the most abundant length inside [size_lo, size_hi] is used.
    """
    n, keys, seqs = _accumulate(fh, chrom, start, end, xm_as_nh=xm_as_nh)
    return _build(n, keys, seqs, start, end, dominant=dominant, size_lo=size_lo,
                  size_hi=size_hi, max_bins=max_bins, top_seqs=top_seqs)


def _empty(start, end):
    return _build(0, Counter(), Counter(), start, end)


_CHUNK = 50  # clusters per worker task


def _chunk_signals(bam_path, items, xm_as_nh):
    """(signals by name, failed count, first error) for plain-tuple clusters.

    items: (name, chrom, start, end, dicer_call). Opens its own AlignmentFile,
    so it can run in a worker process.
    """
    out = {}
    failed, first = 0, None
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        refs = set(fh.references)
        for name, chrom, start, end, dom in items:
            if start is None or end is None:
                failed += 1
                first = first or f"{name}: no coordinates"
                out[name] = dict(_empty(0, 0), failed=True)
                continue
            if chrom not in refs:
                out[name] = _empty(start, end)
                continue
            try:
                n, keys, seqs = _accumulate(fh, chrom, start, end, xm_as_nh=xm_as_nh)
                out[name] = _build(n, keys, seqs, start, end, dominant=dom)
            except Exception as e:
                failed += 1
                first = first or f"{type(e).__name__}: {e}"
                out[name] = dict(_empty(start, end), failed=True)
    return out, failed, first


def signals_for(bam_path, clusters, xm_as_nh=None, threads=1):
    """One figure-data entry per cluster name; a failing region is zero-filled.

    threads > 1 splits the clusters into chunks run in worker processes (each
    opens its own BAM handle); the result is the same dict as threads=1.
    """
    if xm_as_nh is None:
        xm_as_nh = samio.prefers_xm_as_nh(bam_path)
    items = [(c.name, c.chrom, c.start, c.end, getattr(c, "dicer_call", None))
             for c in clusters]
    chunks = [items[i:i + _CHUNK] for i in range(0, len(items), _CHUNK)]
    if threads > 1 and len(chunks) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=min(threads, len(chunks))) as ex:
            results = list(ex.map(_chunk_signals, [bam_path] * len(chunks), chunks,
                                  [xm_as_nh] * len(chunks)))
    else:
        results = [_chunk_signals(bam_path, ch, xm_as_nh) for ch in chunks]
    out, failed, first = {}, 0, None
    for o, f, e in results:  # chunk order keeps the first error deterministic
        out.update(o)
        failed += f
        first = first or e
    if failed:
        utils.log(f"report figures: {failed} cluster(s) failed to read from the BAM "
                  f"(zero-filled); first error: {first}", "WARN")
    return out
