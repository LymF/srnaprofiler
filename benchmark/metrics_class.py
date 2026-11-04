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


def match_loci(loci: list[dict], regions: list[Region],
               min_overlap_frac: float) -> list[tuple[dict, Region | None]]:
    """Pair each locus with its best-overlapping region (or None), in input order.

    A region matches when overlap/locus_len or overlap/region_len >= min_overlap_frac.
    Ties on overlap go to the region that comes first in `regions`.
    """
    by_chrom: dict[str, list[tuple[int, int, int]]] = defaultdict(list)
    for idx, r in enumerate(regions):
        by_chrom[r.chrom].append((r.start, idx, r.end))
    starts: dict[str, list[int]] = {}
    max_len: dict[str, int] = {}
    for chrom, items in by_chrom.items():
        items.sort()
        starts[chrom] = [s for s, _, _ in items]
        max_len[chrom] = max(e - s for s, _, e in items)

    out = []
    for loc in loci:
        chrom, ls, le = loc["chrom"], loc["start"], loc["end"]
        best, best_ov, best_idx = None, 0, -1
        if chrom in by_chrom:
            items = by_chrom[chrom]
            hi = bisect_left(starts[chrom], le)  # regions with start < le
            lo = bisect_left(starts[chrom], ls - max_len[chrom])
            for s, idx, e in items[lo:hi]:
                ov = min(e, le) - max(s, ls)
                if ov <= 0:
                    continue
                if ov / (le - ls) >= min_overlap_frac or ov / (e - s) >= min_overlap_frac:
                    if ov > best_ov or (ov == best_ov and idx < best_idx):
                        best, best_ov, best_idx = regions[idx], ov, idx
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


def classification_metrics(loci: list[dict], regions: list[Region], cfg: Config) -> dict:
    matched = match_loci(loci, regions, cfg.min_overlap_frac)
    scored = [(loc, reg) for loc, reg in matched if loc["cls"] != AMBIGUOUS]

    truth_count: dict[str, int] = defaultdict(int)
    for r in regions:
        truth_count[r.cls] += 1
    classes = set(truth_count) | {loc["cls"] for loc, _ in scored}

    per_class = {}
    for c in sorted(classes):
        tp = fp = 0
        hit_regions = set()
        for loc, reg in scored:
            if loc["cls"] != c or reg is None:
                continue
            if reg.cls == c:
                tp += 1
                hit_regions.add(id(reg))
            else:
                fp += 1
        precision = tp / (tp + fp) if tp + fp else None
        recall = len(hit_regions) / truth_count[c] if truth_count[c] else None
        per_class[c] = {
            "precision": precision,
            "fdr": None if precision is None else 1.0 - precision,
            "recall": recall,
            "tp": tp,
            "fp": fp,
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
