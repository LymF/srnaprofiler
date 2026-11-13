"""Per-locus genomic window: annotation features sliced into compact tracks.

Pure functions, no I/O. Coordinates are 0-based half-open genomic positions;
nothing is clipped to the window (the renderer clips).

Gene track rules: one entry per gene, drawn from its LONGEST transcript (ties
broken by transcript id); every "g" entry has at least one block (a gene with
no transcripts or no blocks gets one "E" block over its span). Exon-like
groups with no transcript row (GTF without transcript lines) attach to the
smallest gene row of the same strand that contains them, else stay their own
entry. The entry's s/e/hl come from the drawn transcript span, so a gene whose
other isoform overlaps the window but whose longest transcript does not is
not listed.
"""
from annotation import LAYER_KINDS, _gene_class
from utils import log

FEATURE_CAP = 60  # entries kept per track; the rest are counted in "mo"

_BLOCK_TYPE = {"exon": "E", "noncoding_exon": "E", "pseudogenic_exon": "E",
               "cds": "C", "five_prime_utr": "U", "three_prime_utr": "U", "utr": "U"}
_FLAT_TRACK = {"te": "t", "ncrna": "r", "phas": "p"}


def window_for(start, end, flank, contig_len):
    """The locus plus flank on each side, clamped to [0, contig_len]."""
    return max(0, start - flank), min(contig_len, end + flank)


def _overlaps(s, e, start, end):
    return s < end and e > start


def _cap(rows, locus, s_i, e_i, hl_i):
    """Sort by (start, end, label), keep highlighted first up to the cap."""
    rows = sorted(rows, key=lambda r: (r[s_i], r[e_i], r[0]))
    if len(rows) <= FEATURE_CAP:
        return rows, 0
    hl = [r for r in rows if r[hl_i]]
    rest = [r for r in rows if not r[hl_i]]
    kept = (hl + rest)[:FEATURE_CAP]
    kept.sort(key=lambda r: (r[s_i], r[e_i], r[0]))
    return kept, len(rows) - FEATURE_CAP


def _fetch(layers, chrom, start, end):
    out = []
    for la in layers:
        idx = la.index.get(chrom)
        if idx is not None:
            out.extend(idx.overlaps(start, end))
    return out


def _gene_entries(feats, ws, we, locus):
    spans = {}   # fid -> span feature
    blocks = []  # (block type, feature)
    anon = []    # span rows without an id
    flat = []    # BED-like rows: one entry, one E block
    for f in feats:
        t = f.ftype.lower()
        if _gene_class(f.ftype) == "span":
            if f.fid:
                spans.setdefault(f.fid, f)
            else:
                anon.append(f)
        elif t == "feature":
            flat.append(f)
        elif t in _BLOCK_TYPE:
            blocks.append((_BLOCK_TYPE[t], f))
        # start/stop codons and anything else: ignored
    by_parent = {}
    for bt, f in blocks:
        by_parent.setdefault(f.parent, []).append((bt, f))
    # a transcript is a span row that carries blocks or hangs under another span row
    tx_ids = {fid for fid, f in spans.items() if fid in by_parent or f.parent in spans}
    for p in by_parent:
        if p and p not in spans:
            tx_ids.add(p)  # blocks under a transcript row that is not in the file
    groups = {}  # gene key -> {"gene": row|None, "tx": {tid: row|None}}
    for fid, f in spans.items():
        if fid in tx_ids:
            key = f.parent or fid
            groups.setdefault(key, {"gene": None, "tx": {}})["tx"][fid] = f
        else:
            groups.setdefault(fid, {"gene": None, "tx": {}})["gene"] = f
    genes = [(key, g["gene"]) for key, g in groups.items() if g["gene"] is not None]
    for tid in sorted(t for t in tx_ids if t not in spans):
        bl = by_parent.get(tid, [])
        gs, ge = min(f.start for _, f in bl), max(f.end for _, f in bl)
        strand = bl[0][1].strand
        home = [(gf.end - gf.start, key) for key, gf in genes
                if gf.strand == strand and gf.start <= gs and gf.end >= ge]
        key = min(home)[1] if home else tid
        groups.setdefault(key, {"gene": None, "tx": {}})["tx"][tid] = None
    for f in anon:
        groups.setdefault(id(f), {"gene": f, "tx": {}})
    rows = []
    for g in groups.values():
        best = None  # (span, tid, row)
        for tid, row in g["tx"].items():
            bl = by_parent.get(tid, [])
            if bl:
                s, e = min(f.start for _, f in bl), max(f.end for _, f in bl)
            elif row is not None:
                s, e = row.start, row.end
            else:
                continue
            if best is None or (e - s, -s) > (best[0][1] - best[0][0], -best[0][0]):
                best = ((s, e), tid, row)
            elif (e - s, -s) == (best[0][1] - best[0][0], -best[0][0]) and tid < best[1]:
                best = ((s, e), tid, row)
        gene = g["gene"]
        if best is None:
            if gene is None:
                continue
            s, e, bls, tr = gene.start, gene.end, [], None
        else:
            (s, e), tid, tr = best
            bls = by_parent.get(tid, [])
        src = gene or tr
        label = (gene.label if gene else "") or (tr.label if tr else "")
        strand = src.strand if src else bls[0][1].strand
        if not _overlaps(s, e, ws, we):
            continue
        bl = sorted(([bt, f.start, f.end] for bt, f in bls),
                    key=lambda b: (b[1], b[2], b[0])) or [["E", s, e]]
        hl = 1 if _overlaps(s, e, locus[0], locus[1]) else 0
        rows.append([label, strand, s, e, hl, bl])
    for bt, f in by_parent.get("", []):  # orphan blocks without any parent
        if _overlaps(f.start, f.end, ws, we):
            hl = 1 if _overlaps(f.start, f.end, locus[0], locus[1]) else 0
            rows.append([f.label, f.strand, f.start, f.end, hl, [[bt, f.start, f.end]]])
    for f in flat:
        if _overlaps(f.start, f.end, ws, we):
            hl = 1 if _overlaps(f.start, f.end, locus[0], locus[1]) else 0
            rows.append([f.label, f.strand, f.start, f.end, hl,
                         [["E", f.start, f.end]]])
    return rows


def window_features(chrom, ws, we, locus, layers):
    """Features of every layer overlapping [ws, we), as compact track rows."""
    out = {"g": [], "t": [], "r": [], "p": [], "mo": {}}
    kinds = {la.kind for la in layers}
    out["ly"] = [k for k in LAYER_KINDS if k in kinds]
    lo, hl_ = locus
    gene_layers = [la for la in layers if la.kind == "gene"]
    feats = _fetch(gene_layers, chrom, ws, we)
    if feats:
        # a gene may extend past the window: refetch its whole span so its
        # transcripts and blocks are complete
        lo_x = min(f.start for f in feats)
        hi_x = max(f.end for f in feats)
        if lo_x < ws or hi_x > we:
            feats = _fetch(gene_layers, chrom, min(ws, lo_x), max(we, hi_x))
    rows = _gene_entries(feats, ws, we, locus)
    out["g"], out["mo"]["g"] = _cap(rows, locus, 2, 3, 4)
    for kind, key in _FLAT_TRACK.items():
        rows = []
        for f in _fetch([la for la in layers if la.kind == kind], chrom, ws, we):
            hl = 1 if _overlaps(f.start, f.end, lo, hl_) else 0
            second = f.name if kind == "te" else f.ftype
            rows.append([f.label, second, f.strand, f.start, f.end, hl])
        out[key], out["mo"][key] = _cap(rows, locus, 3, 4, 5)
    return out


def build_windows(clusters, layers, flank, contig_lengths):
    """{cluster name: window, locus and features} for each cluster passed."""
    res = {}
    warned = set()
    for c in clusters:
        n = contig_lengths.get(c.chrom)
        if n is None:
            n = c.end
            if c.chrom not in warned:
                warned.add(c.chrom)
                log(f"report: contig {c.chrom} has no length; the right flank "
                    f"is clipped to the cluster end", "WARN")
        ws, we = window_for(c.start, c.end, flank, n)
        res[c.name] = {"w": [ws, we], "L": [c.start, c.end],
                       **window_features(c.chrom, ws, we, (c.start, c.end), layers)}
    return res
