"""Classification metrics: match predicted loci to truth regions and score them."""
import csv
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path

import numpy as np

from benchmark.manifest import Config
from benchmark.truth import Region

AMBIGUOUS = "ambiguous"
QUANTILES = {"p5": 5, "p25": 25, "p50": 50, "p75": 75, "p95": 95}
N_BINS = 10


def read_loci(path: Path) -> list[dict]:
    """Parse a clusters.tsv into dicts (start/end int, reads float, confidence float|None)."""
    loci = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["reads"] = float(row["reads"])
            conf = row.get("confidence")
            row["confidence"] = None if conf in (None, "", "NA") else float(conf)
            loci.append(row)
    return loci


class IntervalIndex:
    """Per-chrom index over (chrom, start, end) intervals, sorted by start, for overlap queries."""

    def __init__(self, intervals: list[tuple[str, int, int]]):
        by_chrom: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
        for idx, (chrom, start, end) in enumerate(intervals):
            by_chrom[chrom].append((start, idx, end))
        self._items = by_chrom
        self._starts: dict[str, list[int]] = {}
        self._max_len: dict[str, int] = {}
        for chrom, items in by_chrom.items():
            items.sort()
            self._starts[chrom] = [s for s, _, _ in items]
            self._max_len[chrom] = max(e - s for s, _, e in items)

    def overlapping(self, chrom: str, start: int, end: int):
        """Yield (idx, overlap_bp) for each indexed interval with > 0 bp overlap, by start order."""
        if chrom not in self._items:
            return
        starts = self._starts[chrom]
        hi = bisect_left(starts, end)  # intervals with start < end
        lo = bisect_left(starts, start - self._max_len[chrom])
        for s, idx, e in self._items[chrom][lo:hi]:
            ov = min(e, end) - max(s, start)
            if ov > 0:
                yield idx, ov


def match_loci(loci: list[dict], regions: list[Region],
               min_overlap_frac: float) -> list[tuple[dict, Region | None]]:
    """Pair each locus with its best-overlapping region (or None), in input order.

    A region matches when overlap/locus_len or overlap/region_len >= min_overlap_frac.
    Ties on overlap go to the smaller region, then to the one that comes first in `regions`
    (so a miRNA precursor inside a larger cluster does not depend on the manifest order).
    """
    index = IntervalIndex([(r.chrom, r.start, r.end) for r in regions])
    out = []
    for loc in loci:
        ls, le = loc["start"], loc["end"]
        best, best_ov, best_idx = None, 0, -1
        best_len = 0
        for idx, ov in index.overlapping(loc["chrom"], ls, le):
            r = regions[idx]
            if ov / (le - ls) >= min_overlap_frac or ov / (r.end - r.start) >= min_overlap_frac:
                length = r.end - r.start
                if (best is None or ov > best_ov
                        or (ov == best_ov and (length, idx) < (best_len, best_idx))):
                    best, best_ov, best_idx, best_len = r, ov, idx, length
        out.append((loc, best))
    return out


def _quantiles(values: list[float]) -> dict[str, float]:
    if not values:
        return {k: 0.0 for k in QUANTILES}
    return {k: float(np.percentile(values, q)) for k, q in QUANTILES.items()}


def _ece(pairs: list[tuple[float, bool]]) -> float:
    if not pairs:
        return 0.0
    bins: list[list[tuple[float, bool]]] = [[] for _ in range(N_BINS)]
    for conf, correct in pairs:
        bins[min(int(conf * N_BINS), N_BINS - 1)].append((conf, correct))
    n = len(pairs)
    total = 0.0
    for b in bins:
        if b:
            acc = sum(c for _, c in b) / len(b)
            mean_conf = sum(p for p, _ in b) / len(b)
            total += len(b) / n * abs(acc - mean_conf)
    return total


def _recalled_regions(loci: list[dict], regions: list[Region],
                      min_overlap_frac: float) -> set[int]:
    """id() of each region overlapped by at least one locus (match_loci's overlap criterion)."""
    if not loci or not regions:
        return set()
    index = IntervalIndex([(r.chrom, r.start, r.end) for r in regions])
    hit = set()
    for loc in loci:
        ls, le = loc["start"], loc["end"]
        for idx, ov in index.overlapping(loc["chrom"], ls, le):
            r = regions[idx]
            if ov / (le - ls) >= min_overlap_frac or ov / (r.end - r.start) >= min_overlap_frac:
                hit.add(id(r))
    return hit


def classification_metrics(loci: list[dict], regions: list[Region], cfg: Config,
                           expressed: set[int] | None = None) -> dict:
    """Per-class precision/recall plus unknown and calibration.

    precision is over annotated loci only: loci predicted C that match a truth region
    (tp + fp = n_matched); the n_unknown loci predicted C with no match are outside it.
    Precision is per locus (best-matching region). recall is per truth region: a region of
    class C is recalled when ANY non-ambiguous locus called C overlaps it (same criterion as
    match_loci), so one locus covering several adjacent regions recalls all of them.
    expressed: id() of the truth regions with enough reads (metrics_expr); None = unknown.
    """
    matched = match_loci(loci, regions, cfg.min_overlap_frac)
    scored = [(loc, reg) for loc, reg in matched if loc["cls"] != AMBIGUOUS]

    truth_count: dict[str, int] = defaultdict(int)
    for r in regions:
        truth_count[r.cls] += 1
    classes = set(truth_count) | {loc["cls"] for loc, _ in scored}

    per_class = {}
    for c in sorted(classes):
        tp = fp = n_called = 0
        for loc, reg in scored:
            if loc["cls"] != c:
                continue
            n_called += 1
            if reg is None:
                continue
            if reg.cls == c:
                tp += 1
            else:
                fp += 1
        precision = tp / (tp + fp) if tp + fp else None
        hit_regions = _recalled_regions(
            [loc for loc, _ in scored if loc["cls"] == c],
            [r for r in regions if r.cls == c], cfg.min_overlap_frac)
        recall = len(hit_regions) / truth_count[c] if truth_count[c] else None
        n_expr = n_expr_hit = None
        if expressed is not None:
            n_expr = sum(1 for r in regions if r.cls == c and id(r) in expressed)
            n_expr_hit = len(hit_regions & expressed)
        per_class[c] = {
            "precision": precision,
            "fdr": None if precision is None else 1.0 - precision,
            "recall": recall,
            "recall_expressed": n_expr_hit / n_expr if n_expr else None,
            "n_truth": truth_count[c],
            "n_truth_expressed": n_expr,
            "tp": tp,
            "fp": fp,
            "n_called": n_called,
            "n_matched": tp + fp,
            "n_unknown": n_called - tp - fp,
        }

    unknown = [loc for loc, reg in matched if reg is None]
    ece_pairs = [(loc["confidence"], loc["cls"] == reg.cls) for loc, reg in scored
                 if reg is not None and loc["confidence"] is not None]
    return {
        "per_class": per_class,
        "ambiguous_coverage": len(scored) / len(loci) if loci else 0.0,
        "unknown": {
            "count": len(unknown),
            "length": _quantiles([u["end"] - u["start"] for u in unknown]),
            "reads": _quantiles([u["reads"] for u in unknown]),
        },
        "ece": _ece(ece_pairs),
    }
