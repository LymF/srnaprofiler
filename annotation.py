"""Annotation file parsers (GFF3/GTF/BED/RepeatMasker) and an interval index.

Coordinates are 0-based half-open. Malformed lines are skipped, never raised.
"""
import gzip
import re
from bisect import bisect_left
from collections import namedtuple
from itertools import accumulate

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
