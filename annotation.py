"""Annotation file parsers (GFF3/GTF/BED/RepeatMasker) and an interval index.

Coordinates are 0-based half-open. Malformed lines are skipped, never raised.
"""
import gzip
import re
from bisect import bisect_left
from collections import namedtuple
from itertools import accumulate

from utils import log

Feature = namedtuple("Feature", "start end strand ftype name")

_STRANDS = ("+", "-", ".")
_CLASS_KEYS = ("classification", "Class", "class", "family", "Family")


def _open(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def _strand(s):
    return s if s in _STRANDS else "."


def _attrs_gff3(text):
    out = {}
    for part in text.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def _attrs_gtf(text):
    out = {}
    for m in re.finditer(r'(\w+)\s+(?:"([^"]*)"|([^;\s]+))\s*;?', text):
        out[m.group(1)] = m.group(2) if m.group(2) is not None else m.group(3)
    return out


def _gff_lines(lines, gtf):
    for line in lines:
        line = line.rstrip("\r\n")
        if line.startswith("##FASTA"):
            return
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 8:
            continue
        try:
            start, end = int(cols[3]), int(cols[4])
        except ValueError:
            continue
        if start < 1 or start > end:
            continue
        ftype = cols[2]
        attrs = {}
        if len(cols) > 8:
            attrs = _attrs_gtf(cols[8]) if gtf else _attrs_gff3(cols[8])
        name = ftype
        for key in _CLASS_KEYS:
            if attrs.get(key):
                name = attrs[key]
                break
        yield cols[0], Feature(start - 1, end, _strand(cols[6]), ftype, name)


def _bed_lines(lines):
    for line in lines:
        line = line.rstrip("\r\n")
        if (not line.strip() or line.startswith(("#", "track", "browser"))):
            continue
        cols = line.split("\t") if "\t" in line else line.split()
        if len(cols) < 3:
            continue
        try:
            start, end = int(cols[1]), int(cols[2])
        except ValueError:
            continue
        if start < 0 or start > end:
            continue
        strand = _strand(cols[5]) if len(cols) > 5 else "."
        yield cols[0], Feature(start, end, strand, "feature", "")


def _rm_lines(lines):
    for line in lines:
        cols = line.split()
        if len(cols) < 11:
            continue
        try:
            begin, end = int(cols[5]), int(cols[6])
        except ValueError:
            # header lines ("begin", "end") and blanks land here
            continue
        if begin < 1 or begin > end:
            continue
        strand = "-" if cols[8] == "C" else "+" if cols[8] == "+" else "."
        yield cols[4], Feature(begin - 1, end, strand, "feature", cols[10])


def parse_annotation(path):
    """Parse an annotation file into {contig: [Feature, ...]}."""
    lower = str(path).lower()
    if lower.endswith(".gz"):
        lower = lower[:-3]
    out = {}
    with _open(path) as fh:
        if lower.endswith(".out"):
            gen = _rm_lines(fh)
        elif lower.endswith(".bed"):
            gen = _bed_lines(fh)
        elif lower.endswith(".gtf"):
            gen = _gff_lines(fh, gtf=True)
        else:
            gen = _gff_lines(fh, gtf=False)
        for contig, feat in gen:
            out.setdefault(contig, []).append(feat)
    return out


class IntervalIndex:
    """Overlap queries over one contig's features."""

    def __init__(self, features):
        self._feats = sorted(features, key=lambda f: (f.start, f.end))
        self._starts = [f.start for f in self._feats]
        self._maxend = list(accumulate((f.end for f in self._feats), max))

    def overlaps(self, start, end):
        """Features overlapping [start, end)."""
        hi = bisect_left(self._starts, end)  # features with start < end
        res = []
        i = hi - 1
        while i >= 0 and self._maxend[i] > start:
            if self._feats[i].end > start:
                res.append(self._feats[i])
            i -= 1
        res.reverse()
        return res


LAYER_KINDS = ("gene", "te", "ncrna", "phas")
_SPAN_TYPES = {"gene", "mRNA", "transcript", "ncRNA_gene"}


class Layer:
    """One annotation layer: a kind label plus per-contig interval indexes."""

    def __init__(self, kind, path, index, n_features):
        self.kind = kind
        self.path = path
        self.index = index
        self.n_features = n_features


def parse_layer_spec(spec):
    """Split 'kind:path' on the first ':' and validate the kind."""
    kind, sep, path = spec.partition(":")
    if kind not in LAYER_KINDS or not sep or not path:
        raise ValueError(
            f"invalid annotation layer {spec!r}: expected <kind>:<path> with "
            f"kind one of {', '.join(LAYER_KINDS)}")
    return kind, path


def load_layers(specs, known_contigs):
    """Load layers, dropping contigs absent from known_contigs (warned once per layer)."""
    layers = []
    for spec in specs:
        kind, path = parse_layer_spec(spec)
        parsed = parse_annotation(path)
        index = {}
        n_kept = 0
        skipped = 0
        unmatched = []
        for contig, feats in parsed.items():
            if contig in known_contigs:
                index[contig] = IntervalIndex(feats)
                n_kept += len(feats)
            else:
                skipped += len(feats)
                unmatched.append(contig)
        if skipped:
            log(f"annotation layer {kind} ({path}): skipped {skipped} "
                f"features on {len(unmatched)} contigs not in the reference "
                f"(e.g. {', '.join(unmatched[:3])})", "WARN")
        layers.append(Layer(kind, path, index, n_kept))
    return layers


def _covered(feats, start, end):
    """Bases of [start, end) covered by the union of feats."""
    ivs = sorted((max(f.start, start), min(f.end, end)) for f in feats)
    total = 0
    cur_s = cur_e = None
    for a, b in ivs:
        if b <= a:
            continue
        if cur_e is None or a > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = a, b
        else:
            cur_e = max(cur_e, b)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


def _relation(locus_strand, feat_strand):
    if locus_strand not in ("+", "-"):
        return "."
    if feat_strand not in ("+", "-"):
        return "both"
    return "sense" if locus_strand == feat_strand else "antisense"


def _best(feats, start, end):
    return max(feats, key=lambda f: min(f.end, end) - max(f.start, start))


def _entry(kind, feats, start, end, strand, frac):
    b = _best(feats, start, end)
    return f"{kind}:{b.name or b.ftype}:{_relation(strand, b.strand)}:{frac:.2f}"


def summarize_locus(chrom, start, end, strand, layers, min_overlap):
    """Contexts and a text summary of what a locus overlaps, per layer."""
    length = end - start
    contexts = set()
    entries = []
    annotated = False  # contig present in some gene layer
    for layer in layers:
        if layer.kind == "gene" and chrom in layer.index:
            annotated = True
        idx = layer.index.get(chrom)
        if idx is None or length <= 0:
            continue
        feats = idx.overlaps(start, end)
        if not feats:
            continue
        if layer.kind == "gene":
            exon_like = [f for f in feats if f.ftype not in _SPAN_TYPES]
            span = [f for f in feats if f.ftype in _SPAN_TYPES]
            ef = _covered(exon_like, start, end) / length
            sf = _covered(span, start, end) / length
            if ef >= min_overlap:
                contexts.add("gene_exon")
                entries.append(_entry("gene", exon_like, start, end, strand, ef))
            elif sf >= min_overlap:
                contexts.add("gene_intron")
                entries.append(_entry("gene", span, start, end, strand, sf))
        else:
            frac = _covered(feats, start, end) / length
            if frac >= min_overlap:
                contexts.add(layer.kind)
                entries.append(_entry(layer.kind, feats, start, end, strand, frac))
    if annotated and not contexts:
        contexts.add("intergenic")
    return {"contexts": contexts, "text": "|".join(entries)}


def support_for(cls, contexts, table):
    """'supports', 'conflicts' or 'neutral' for a class given a locus's contexts."""
    spec = table.get(cls)
    if not spec or not contexts:
        return "neutral"
    if contexts & set(spec.get("supports") or ()):
        return "supports"
    if contexts & set(spec.get("conflicts") or ()):
        return "conflicts"
    return "neutral"


def apply_annotation(records, layers, cfg):
    """Fill annotation fields; nudge confidence and flag review. Never changes cls."""
    counts = {"none": 0, "supports": 0, "conflicts": 0, "neutral": 0}
    if not layers:
        counts["none"] = len(records)
        return counts
    nudge = cfg.get("nudge", 0.05)
    min_overlap = cfg.get("min_overlap", 0.5)
    table = cfg.get("classes") or {}
    for r in records:
        res = summarize_locus(r.chrom, r.start, r.end, r.strand_call, layers, min_overlap)
        r.annotation = res["text"]
        sup = support_for(r.cls, res["contexts"], table)
        r.annotation_support = sup
        counts[sup] += 1
        if r.cls in table and r.cls not in ("miRNA", "ambiguous"):
            if sup == "supports":
                r.confidence = min(1.0, round(r.confidence + nudge, 3))
            elif sup == "conflicts":
                r.confidence = max(0.0, round(r.confidence - nudge, 3))
                r.review = True
    return counts
