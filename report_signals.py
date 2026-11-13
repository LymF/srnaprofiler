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
READ_ROWS = 25  # read variants drawn per hairpin
_SKIP = 0x4 | 0x100 | 0x800  # unmapped, secondary, supplementary (samio._keep)


def _r(x):
    return round(float(x), 3)


class HairpinReads:
    """Reads inside one miRNA precursor, collected while the cluster region is read.

    spec: (lo, hi, strand, precursor RNA sequence, mature (a, b) | None, star | None); lo, hi
    are the genomic half-open window of the precursor, a and b positions in the precursor
    (strand-oriented: on "-" position 0 is the highest genomic coordinate). Only reads of the
    precursor strand that lie fully inside the window and align without indels or clips are
    kept, collapsed by (a, b, sequence) with 1/NH-weighted counts. The sequence is written in
    the precursor's orientation as RNA, lower case where it matches the precursor and upper
    case at a mismatch.
    """

    def __init__(self, spec):
        self.lo, self.hi, self.strand, self.prec, self.m, self.s = spec
        self.rows = Counter()
        self.anti = self.ovh = 0.0

    def add(self, rec, nh):
        if rec.reference_end <= self.lo or rec.reference_start >= self.hi:
            return   # elsewhere in the cluster, not part of this precursor
        w = 1.0 / nh
        if rec.is_reverse != (self.strand == "-"):
            self.anti += w
            return
        if rec.reference_start < self.lo or rec.reference_end > self.hi:
            self.ovh += w
            return
        if any(op not in (0, 7, 8) for op, _ in rec.cigartuples or ()):
            self.ovh += w
            return
        seq = rec.query_sequence.upper()
        if rec.is_reverse:
            seq = revcomp(seq)
        seq = seq.replace("T", "U")
        a, b = (self.hi - rec.reference_end, self.hi - rec.reference_start) \
            if self.strand == "-" else (rec.reference_start - self.lo, rec.reference_end - self.lo)
        ref = self.prec[a:b]
        self.rows[a, b, "".join(x.lower() if x == y else x for x, y in zip(seq, ref))] += w

    def _kind(self, a, b):
        for key, seg in (("m", self.m), ("s", self.s)):
            if seg and min(b, seg[1]) - max(a, seg[0]) >= 0.5 * (b - a):
                return key
        return ""

    def build(self, top=READ_ROWS):
        """{"rd": top rows [a, b, sequence, reads, kind], "cov": per-base depth, "more": distinct
        sequences not in rd, "mr": their reads, "anti", "ovh": reads left out}."""
        L = len(self.prec)
        diff = [0.0] * (L + 1)
        for (a, b, _), w in self.rows.items():
            diff[a] += w
            diff[b] -= w
        cov, run = [], 0.0
        for i in range(L):
            run += diff[i]
            cov.append(_r(run))
        # ties at 1e-9 so float summation order cannot reorder them
        order = sorted(self.rows.items(), key=lambda kv: (-round(kv[1], 9), kv[0]))
        return {"rd": [[a, b, q, _r(w), self._kind(a, b)] for (a, b, q), w in order[:top]],
                "cov": cov, "more": max(0, len(order) - top),
                "mr": _r(sum(w for _, w in order[top:])),
                "anti": _r(self.anti), "ovh": _r(self.ovh)}


def _accumulate(fh, chrom, start, end, *, xm_as_nh, hp=None):
    """One pass over the region: (records kept, key Counter, sequence Counter).

    keys: (strand, pos, length, nh, 5' base) -> record count.
    seqs: (is_reverse, stored sequence) -> 1/NH-weighted reads.
    hp: a HairpinReads fed in the same pass; the region read is then the union of the
    cluster window and the precursor window, and the counters still see only the reads
    that overlap the cluster window.
    """
    keys, seqs = Counter(), Counter()
    n = 0
    if hp is None:
        lo, hi = start, end
    else:
        lo, hi = min(start, hp.lo), max(end, hp.hi)
    for rec in fh.fetch(chrom, lo, hi):
        if rec.flag & _SKIP:
            continue
        seq = rec.query_sequence
        if not seq:
            continue
        nh = samio._nh(rec, xm_as_nh)
        if hp is not None:
            hp.add(rec, nh)
            if rec.reference_start >= end or rec.reference_end <= start:
                continue
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

    items: (name, chrom, start, end, dicer_call[, hairpin spec]). Opens its own
    AlignmentFile, so it can run in a worker process. A cluster with a hairpin spec gets
    "hp" (HairpinReads.build) in its entry, read in the same pass as its figure data.
    """
    out = {}
    failed, first = 0, None
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        refs = set(fh.references)
        for name, chrom, start, end, dom, *rest in items:
            if start is None or end is None:
                failed += 1
                first = first or f"{name}: no coordinates"
                out[name] = dict(_empty(0, 0), failed=True)
                continue
            if chrom not in refs:
                out[name] = _empty(start, end)
                continue
            try:
                hp = HairpinReads(rest[0]) if rest and rest[0] else None
                n, keys, seqs = _accumulate(fh, chrom, start, end, xm_as_nh=xm_as_nh, hp=hp)
                out[name] = _build(n, keys, seqs, start, end, dominant=dom)
                if hp is not None:
                    out[name]["hp"] = hp.build()
            except Exception as e:
                failed += 1
                first = first or f"{type(e).__name__}: {e}"
                out[name] = dict(_empty(start, end), failed=True)
    return out, failed, first


def signals_for(bam_path, clusters, xm_as_nh=None, threads=1, hairpins=None):
    """One figure-data entry per cluster name; a failing region is zero-filled.

    threads > 1 splits the clusters into chunks run in worker processes (each
    opens its own BAM handle); the result is the same dict as threads=1.
    hairpins: {cluster name: hairpin spec} for the miRNA clusters (see HairpinReads).
    """
    if xm_as_nh is None:
        xm_as_nh = samio.prefers_xm_as_nh(bam_path)
    hairpins = hairpins or {}
    items = [(c.name, c.chrom, c.start, c.end, getattr(c, "dicer_call", None),
              hairpins.get(c.name)) for c in clusters]
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


# ---- whole-library read lengths (one extra pass, by contig) ----

LEN_LO, LEN_HI = 16, 40  # exact bins; below and above go to the two overflow bins


def _len_slot(n):
    return 0 if n < LEN_LO else LEN_HI - LEN_LO + 2 if n > LEN_HI else n - LEN_LO + 1


def _contigs_lengths(bam_path, contigs, xm_as_nh):
    """1/NH-weighted read counts per length slot over the given contigs."""
    counts = [0.0] * (LEN_HI - LEN_LO + 3)
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        for contig in contigs:
            for rec in fh.fetch(contig):
                if rec.flag & _SKIP:
                    continue
                n = rec.query_length
                if n:
                    counts[_len_slot(n)] += 1.0 / samio._nh(rec, xm_as_nh)
    return counts


def _balanced(stats, parts):
    """Split (contig, alignments) pairs into at most `parts` groups of similar size."""
    groups = [[] for _ in range(max(1, min(parts, len(stats))))]
    load = [0] * len(groups)
    for contig, n in sorted(stats, key=lambda x: -x[1]):
        i = load.index(min(load))
        groups[i].append(contig)
        load[i] += n
    return [g for g in groups if g]


def library_lengths(bam_path, xm_as_nh=None, threads=1):
    """Read-length histogram of the whole library.

    {"lo", "hi", "counts": [<lo, lo..hi, >hi] (len 27), "unmapped": int}. Primary mapped
    alignments, each weighted by 1/NH like the cluster counts (so the total is on the same
    scale as the library size behind RPM); secondary, supplementary and unmapped records are
    not counted. `unmapped` is the BAM's own unmapped record count (from the index), not
    weighted. threads > 1 splits the contigs over worker processes.
    """
    if xm_as_nh is None:
        xm_as_nh = samio.prefers_xm_as_nh(bam_path)
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        stats = fh.get_index_statistics()
        unmapped = fh.nocoordinate + sum(s.unmapped for s in stats)
        work = [(s.contig, s.mapped) for s in stats if s.mapped > 0]
    groups = _balanced(work, threads)
    if threads > 1 and len(groups) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=len(groups)) as ex:
            parts = list(ex.map(_contigs_lengths, [bam_path] * len(groups), groups,
                                [xm_as_nh] * len(groups)))
    else:
        parts = [_contigs_lengths(bam_path, g, xm_as_nh) for g in groups]
    counts = [sum(col) for col in zip(*parts)] if parts else [0.0] * (LEN_HI - LEN_LO + 3)
    return {"lo": LEN_LO, "hi": LEN_HI, "counts": counts, "unmapped": int(unmapped)}
