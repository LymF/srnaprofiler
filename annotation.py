"""Annotation file parsers (GFF3/GTF/BED/RepeatMasker) and an interval index.

Coordinates are 0-based half-open. Malformed lines are skipped, never raised.
"""
import gzip
import os
import re
import zlib
from sys import intern
from urllib.parse import unquote
from bisect import bisect_left
from collections import namedtuple

from utils import log

Feature = namedtuple("Feature", "start end strand ftype name fid parent label",
                    defaults=("", "", ""))

_STRANDS = ("+", "-", ".")
_ID_KEYS_GFF3 = ("ID",)
_ID_KEYS_GTF = ("transcript_id", "gene_id")
_LABEL_KEYS = ("Name", "gene_name", "gene", "locus_tag")
_GTF_GENE_TYPES = ("gene",)
_GTF_TRANSCRIPT_TYPES = ("transcript", "mRNA", "ncRNA", "lnc_RNA", "tRNA",
                         "rRNA", "miRNA", "snRNA", "snoRNA", "primary_transcript")
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
            out[k.strip()] = unquote(v.strip())
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
        fid = ""
        for key in (_ID_KEYS_GTF if gtf else _ID_KEYS_GFF3):
            if attrs.get(key):
                fid = attrs[key]
                break
        if gtf:
            if ftype in _GTF_GENE_TYPES:
                parent = ""
            elif ftype in _GTF_TRANSCRIPT_TYPES:
                parent = attrs.get("gene_id", "")
            else:
                parent = attrs.get("transcript_id", "")
        else:
            parent = attrs.get("Parent", "").split(",")[0].strip()
        label = ""
        for key in _LABEL_KEYS + (("gene_id",) if gtf else ("ID",)):
            if attrs.get(key):
                label = attrs[key]
                break
        yield cols[0], Feature(start - 1, end, _strand(cols[6]), intern(ftype),
                               intern(name), fid, intern(parent), intern(label))


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
        name = cols[3] if len(cols) > 3 and cols[3] != "." else ""
        yield cols[0], Feature(start, end, strand, "feature", name, "", "", intern(name))


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
        yield cols[4], Feature(begin - 1, end, strand, "feature", cols[10], "", "",
                               intern(cols[9]))


def parse_annotation(path, known_contigs=None, unknown=None):
    """Parse an annotation file into {contig: [Feature, ...]}.

    With known_contigs, features on other contigs are not retained; they are
    only counted into the `unknown` dict ({contig: n}) when one is given.
    """
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
            if known_contigs is not None and contig not in known_contigs:
                if unknown is not None:
                    unknown[contig] = unknown.get(contig, 0) + 1
                continue
            out.setdefault(contig, []).append(feat)
    return out


LONG_FEATURE = 100_000  # bp; longer features are scanned linearly, not bisected


class IntervalIndex:
    """Overlap queries over one contig's features.

    Features longer than LONG_FEATURE (chromosome-long records) live in a
    separate list scanned on every query; the rest are bisected, and the
    backward walk stops at start - (longest short feature).
    """

    def __init__(self, features):
        feats = sorted(features, key=lambda f: (f.start, f.end))
        self._long = [f for f in feats if f.end - f.start > LONG_FEATURE]
        self._feats = [f for f in feats if f.end - f.start <= LONG_FEATURE]
        self._starts = [f.start for f in self._feats]
        self._max_len = max((f.end - f.start for f in self._feats), default=0)

    def overlaps(self, start, end):
        """Features overlapping [start, end), ordered by (start, end)."""
        hi = bisect_left(self._starts, end)  # features with start < end
        lo = bisect_left(self._starts, start - self._max_len)
        res = [f for f in self._long if f.start < end and f.end > start]
        res.extend(f for f in self._feats[lo:hi] if f.end > start)
        if self._long:
            res.sort(key=lambda f: (f.start, f.end))
        return res


LAYER_KINDS = ("gene", "te", "ncrna", "phas")
_EXON_TYPES = {"exon", "cds", "five_prime_utr", "three_prime_utr", "utr", "start_codon",
               "stop_codon", "noncoding_exon", "pseudogenic_exon", "feature"}
_SPAN_EXACT = {"gene", "pseudogene", "ncrna_gene", "transcript", "primary_transcript"}


def _gene_class(ftype):
    """'exon', 'span' or None (ignored) for a feature type in a gene layer.

    A whitelist: stock files carry region/chromosome/match/transposable_element
    rows that must never count as exons. Any type containing "transposable_element"
    (also transposable_element_gene, FlyBase) is ignored even though *_gene would
    otherwise be a span. BED rows (ftype "feature") are exon-like.
    """
    t = ftype.lower()
    if "transposable_element" in t:       # dm6 FlyBase TE insertions: the te layer's business, not a gene
        return None
    if t in _EXON_TYPES:
        return "exon"
    if (t in _SPAN_EXACT or t.endswith("_gene") or t.endswith("rna")
            or "transcript" in t):
        return "span"
    return None


# RepeatMasker also reports repeats that are not transposons. They must not count as TE
# evidence (a locus on a (TTTTA)n microsatellite is not a TE-derived siRNA).
_NON_TE_CLASSES = frozenset({"simple_repeat", "low_complexity", "rrna", "trna", "snrna",
                             "scrna", "srprna", "artefact"})


def _is_non_te_repeat(f):
    """True when a te-layer feature is a simple/low-complexity repeat or a structural RNA."""
    for text in (f.name, f.label):
        if text and text.split("/", 1)[0].strip().lower() in _NON_TE_CLASSES:
            return True
    return False


def _is_te_type(ftype):
    return "transposable_element" in ftype.lower()


def _drop_te_trees(feats):
    """Drop transposable_element* rows of a gene layer AND their descendants.

    Araport11 lists TEs as transposable_element_gene -> mRNA -> exon/CDS; the
    children carry ordinary types, so they are found through the parent ids,
    transitively and independent of file order.
    """
    dead = {f.fid for f in feats if _is_te_type(f.ftype) and f.fid}
    if not dead:
        return feats
    kids = [f for f in feats if f.parent and not _is_te_type(f.ftype)]
    gone = set()
    changed = True
    while changed:
        changed = False
        for i, f in enumerate(kids):
            if i not in gone and f.parent in dead:
                gone.add(i)
                if f.fid and f.fid not in dead:
                    dead.add(f.fid)
                changed = True
    drop_ids = {id(kids[i]) for i in gone}
    return [f for f in feats if not _is_te_type(f.ftype) and id(f) not in drop_ids]


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


def check_layer_specs(specs):
    """Fail on a bad label or a missing file, before any heavy work (shared by classify and report)."""
    for spec in specs:
        _, path = parse_layer_spec(spec)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"annotation file not found: {path}")


def load_layers(specs, known_contigs):
    """Load layers, keeping only contigs in known_contigs (the BAM references).

    Contigs of the file absent from the reference are warned about once per
    layer; when NO contig of the layer matches the reference the warning says
    so loudly instead (usually chr1 vs 1 naming).
    """
    layers = []
    for spec in specs:
        kind, path = parse_layer_spec(spec)
        unknown = {}
        try:
            parsed = parse_annotation(path, known_contigs, unknown)
        except (EOFError, OSError, zlib.error) as exc:  # truncated/fake .gz, I/O
            raise ValueError(
                f"cannot read annotation file {path}: {type(exc).__name__}: {exc}"
            ) from exc
        index = {}
        n_kept = 0
        ignored = 0
        non_te = 0
        for contig, feats in parsed.items():
            if kind == "gene":
                n_all = len(feats)
                feats = _drop_te_trees(feats)
                feats = [f for f in feats if _gene_class(f.ftype)]
                ignored += n_all - len(feats)
            if kind == "te":
                n_all = len(feats)
                feats = [f for f in feats if not _is_non_te_repeat(f)]
                non_te += n_all - len(feats)
            index[contig] = IntervalIndex(feats)
            n_kept += len(feats)
        skipped = sum(unknown.values())
        eg = ", ".join(list(unknown)[:3])
        if not parsed and unknown:
            log(f"annotation layer {kind} ({path}): NO contig of this file "
                f"matches the reference ({len(unknown)} contigs, {skipped} "
                f"features, e.g. {eg}); check contig naming, e.g. chr1 vs 1. "
                f"This layer will contribute nothing", "WARN")
        elif skipped:
            log(f"annotation layer {kind} ({path}): skipped {skipped} "
                f"features on {len(unknown)} contigs not in the reference "
                f"(e.g. {eg})", "WARN")
        elif not parsed:
            log(f"annotation layer {kind} ({path}): no features parsed from "
                f"this file", "WARN")
        if non_te:
            log(f"annotation layer {kind} ({path}): left out {non_te} simple/low-complexity "
                f"repeats and structural-RNA repeats (not transposons)", "INFO")
        if ignored:
            log(f"annotation layer {kind} ({path}): ignored {ignored} features "
                f"of non gene/exon types (region, chromosome, match, transposable_element, ...)", "DEBUG")
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


def _clean(text):
    for ch in ("\t", "\r", "\n"):
        text = text.replace(ch, " ")
    return text.replace("|", "_").replace(":", "_")


def _is_gene_level(f):
    t = f.ftype.lower()
    return t in ("gene", "pseudogene") or t.endswith("_gene")


def _gene_name(feats, b, start, end):
    """Name of a gene-layer entry: the label of the gene-level feature the chosen
    feature belongs to (exon -> transcript -> gene through `parent`). Without a
    usable parent chain (BED, GTF without parents) the gene-level feature with the
    largest overlap decides (ties: smaller start, then label); then the chosen
    feature's own label, name or type."""
    by_id = {}
    for f in feats:
        if f.fid and _gene_class(f.ftype) == "span":
            by_id.setdefault(f.fid, f)
    top, cur, seen = None, b, set()
    while cur is not None:
        if _is_gene_level(cur):
            top = cur
            break
        nxt = by_id.get(cur.parent) if cur.parent and cur.parent not in seen else None
        seen.add(cur.parent)
        if nxt is None:
            top = cur if cur is not b else None   # chain ends at a transcript row
            break
        cur = nxt
    if top is not None:
        return top.label or b.label or b.name or b.ftype
    if not b.parent or b.parent not in by_id:
        genes = [f for f in feats if f.label and _is_gene_level(f)]
        if genes:
            return min(genes, key=lambda f: (-(min(f.end, end) - max(f.start, start)),
                                             f.start, f.label)).label
    return b.label or b.name or b.ftype


def _entry(kind, feats, start, end, strand, frac, names=None):
    """kind:name:relation:fraction[:feature_id]; the id only if it is parseable.

    A gene-layer entry carries the feature's label (Name) when present instead of its
    type (`gene:HA8:...`); the other layers keep the family or class name."""
    b = _best(feats, start, end)
    name = _gene_name(names, b, start, end) if kind == "gene" else b.name or b.ftype
    out = f"{kind}:{_clean(name)}:{_relation(strand, b.strand)}:{frac:.2f}"
    fid = getattr(b, "fid", "")
    if fid and "|" not in fid and ":" not in fid:
        out += f":{fid}"
    return out


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
            exon_like = [f for f in feats if _gene_class(f.ftype) == "exon"]
            span = [f for f in feats if _gene_class(f.ftype) == "span"]
            ef = _covered(exon_like, start, end) / length
            sf = _covered(span, start, end) / length
            if ef >= min_overlap:
                contexts.add("gene_exon")
                entries.append(_entry("gene", exon_like, start, end, strand, ef, feats))
            elif sf >= min_overlap:
                contexts.add("gene_intron")
                entries.append(_entry("gene", span, start, end, strand, sf, feats))
        else:
            frac = _covered(feats, start, end) / length
            if frac >= min_overlap:
                contexts.add(layer.kind)
                entries.append(_entry(layer.kind, feats, start, end, strand, frac))
    if annotated and not contexts:
        contexts.add("intergenic")
        entries.append("gene:intergenic:.:1.00")
    return {"contexts": contexts, "text": "|".join(entries)}


_CONTEXT_LAYER = {"te": "te", "ncrna": "ncrna", "phas": "phas", "gene_exon": "gene",
                  "gene_intron": "gene", "intergenic": "gene"}


_ABSENCE_CONFLICTS = {"gene_exon", "gene_intron", "intergenic"}


def _needed_layers(spec):
    """Layer kinds a `conflicts` verdict depends on.

    Only an absence-type conflict (exon with NO TE) claims the absence of the
    supporting feature, so it needs every layer kind that could have supplied
    `supports`. A presence-type conflict (miRNA overlapping a TE) needs nothing
    beyond the layer that produced the context.
    """
    if not set(spec.get("conflicts") or ()) & _ABSENCE_CONFLICTS:
        return set()
    return {_CONTEXT_LAYER[k] for k in (spec.get("supports") or ()) if k in _CONTEXT_LAYER}


def support_for(cls, contexts, table, loaded=None):
    """'supports', 'conflicts' or 'neutral' for a class given a locus's contexts.

    A conflict claims the absence of the supporting feature, so it only counts
    when every layer kind that could have supplied it was loaded (`loaded`
    is a set of layer kinds; None means all of them).
    """
    spec = table.get(cls)
    if not spec or not contexts:
        return "neutral"
    if contexts & set(spec.get("supports") or ()):
        return "supports"
    if contexts & set(spec.get("conflicts") or ()):
        if loaded is not None and not _needed_layers(spec) <= set(loaded):
            return "neutral"
        return "conflicts"
    return "neutral"


CONTEXT_KINDS = ("te", "gene_exon", "gene_intron", "intergenic", "ncrna", "phas")


def _number(cfg, key, default, ok, rule):
    v = cfg.get(key, default)
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not ok(v):
        raise ValueError(f"annotation.{key} must be {rule}, got {v!r}")
    return v


def validate_annotation_config(cfg):
    """Return the annotation block normalised; ValueError naming the bad key.

    `annotation:` null in a YAML is the empty block. Warns once when the block
    has no classes (e.g. a config copied from an old rules.yaml): the layers
    would be reported but never support or conflict with anything.
    """
    if cfg is None:
        cfg = {}
    if not isinstance(cfg, dict):
        raise ValueError(f"annotation must be a mapping, got {type(cfg).__name__}")
    out = dict(cfg)
    _number(cfg, "nudge", 0.05, lambda v: 0 <= v <= 0.2, "a number in [0, 0.2]")
    _number(cfg, "min_overlap", 0.5, lambda v: 0 < v <= 1, "a number in (0, 1]")
    classes = cfg.get("classes")
    if classes is None:
        classes = {}
    if not isinstance(classes, dict):
        raise ValueError("annotation.classes must be a mapping of class -> "
                         "{supports, conflicts}")
    for cls, spec in classes.items():
        if spec is None:
            continue
        if not isinstance(spec, dict):
            raise ValueError(f"annotation.classes.{cls} must be a mapping with "
                             f"supports/conflicts lists")
        for key in ("supports", "conflicts"):
            v = spec.get(key)
            if v is None:
                continue
            if not isinstance(v, (list, tuple)):
                hint = f" (write [{v}] for one context)" if isinstance(v, str) else ""
                raise ValueError(f"annotation.classes.{cls}.{key} must be a "
                                 f"list, got {v!r}{hint}")
            for kind in v:
                if kind not in CONTEXT_KINDS:
                    raise ValueError(
                        f"annotation.classes.{cls}.{key}: unknown context "
                        f"{kind!r}; use one of {', '.join(CONTEXT_KINDS)}")
    out["classes"] = classes
    if not classes:
        log("annotation block has no classes: annotation is reported per locus "
            "but never supports or conflicts with a class (copy the annotation "
            "block from the shipped rules.yaml)", "WARN")
    return out


def apply_annotation(records, layers, cfg):
    """Fill annotation fields; nudge confidence and flag review. Never changes cls."""
    cfg = cfg or {}
    counts = {"none": 0, "supports": 0, "conflicts": 0, "neutral": 0}
    if not layers:
        counts["none"] = len(records)
        return counts
    nudge = cfg.get("nudge", 0.05)
    min_overlap = cfg.get("min_overlap", 0.5)
    table = cfg.get("classes") or {}
    loaded = {la.kind for la in layers}
    for cls, spec in table.items():
        missing = sorted(_needed_layers(spec) - loaded)
        if spec.get("conflicts") and missing:
            log(f"annotation: {cls}: conflicts disabled, no {', '.join(missing)} "
                f"layer loaded (an absent feature is not evidence)", "INFO")
    for r in records:
        res = summarize_locus(r.chrom, r.start, r.end, r.strand_call, layers, min_overlap)
        r.annotation = res["text"]
        sup = support_for(r.cls, res["contexts"], table, loaded)
        r.annotation_support = sup
        counts[sup] += 1
        if r.cls in table and r.cls not in ("miRNA", "ambiguous"):
            if sup == "supports":
                r.confidence = min(1.0, round(r.confidence + nudge, 3))
            elif sup == "conflicts":
                r.confidence = max(0.0, round(r.confidence - nudge, 3))
                r.review = True
    return counts
