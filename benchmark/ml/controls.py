"""Negative-control membership: does a locus fall on an rRNA / protein-coding exon interval."""
from pathlib import Path

import numpy as np
import pandas as pd

from benchmark.metrics_class import IntervalIndex


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.split("\t")
            out.append((f[0], int(f[1]), int(f[2])))
    return out


def overlaps(index: IntervalIndex, intervals: list[tuple[str, int, int]], chrom: str,
             start: int, end: int, min_frac: float):
    """Yield idx of indexed intervals matching under the match_loci rule (>= min_frac of either)."""
    for idx, ov in index.overlapping(chrom, start, end):
        _, s, e = intervals[idx]
        if ov / (end - start) >= min_frac or ov / (e - s) >= min_frac:
            yield idx


def is_control(table: pd.DataFrame, controls_bed: Path,
               min_overlap_frac: float | None = None) -> pd.Series:
    """True where a locus overlaps a control interval by >= min_overlap_frac of the locus or of
    the interval (the match_loci rule; default from ml/config.yaml)."""
    if min_overlap_frac is None:
        from benchmark.ml.dataset import load_ml_config
        min_overlap_frac = float(load_ml_config()["min_overlap_frac"])
    ctrl = read_bed(Path(controls_bed))
    index = IntervalIndex(ctrl) if ctrl else None
    flags = np.zeros(len(table), dtype=bool)
    if index is not None:
        for i, (c, s, e) in enumerate(zip(table["chrom"], table["start"], table["end"])):
            flags[i] = next(overlaps(index, ctrl, c, int(s), int(e), min_overlap_frac), None) is not None
    return pd.Series(flags, index=table.index, dtype=bool)
