"""Independent biology of unannotated calls (selection criterion v2, condition 3).

A call is *supported* when evidence the gate and the features do not use backs it:
  * the locus recurs in at least half of the other libraries of its (organism, tissue) group
    (`repro_frac` >= 0.5; NaN, e.g. a single library, is not supported);
  * for class piRNA, and only for organisms with a TE file configured (ml/config.yaml `te`),
    it also overlaps a RepeatMasker transposon interval by >= 1 bp.
Organisms without a TE file (Bombyx, plants for now) are judged by reproducibility alone.
Only unannotated calls count: label null, not carved, prediction neither ambiguous nor null.
Everything else is False and outside the denominator (see `unannotated_calls`).
"""
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from benchmark.metrics_class import IntervalIndex
from benchmark.ml.controls import read_bed

log = logging.getLogger(__name__)
REPRO_MIN = 0.5
AMBIGUOUS = "ambiguous"
TRUTH_DIR = Path(__file__).resolve().parents[1] / "truth"


class TEIndex:
    """Transposon intervals (chrom, start0, end) with an overlap query."""

    def __init__(self, intervals):
        self.intervals = list(intervals)
        self.index = IntervalIndex(self.intervals) if self.intervals else None

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        if self.index is None:
            return False
        return next(self.index.overlapping(chrom, int(start), int(end)), None) is not None


def load_te(path: Path) -> TEIndex:
    return TEIndex(read_bed(Path(path)))


def load_te_indexes(cfg: dict, truth_dir: Path = TRUTH_DIR) -> dict[str, TEIndex]:
    """organism -> TEIndex for the organisms whose configured TE bed exists."""
    out = {}
    for org, rel in (cfg.get("te") or {}).items():
        path = Path(truth_dir) / rel
        if path.exists():
            out[org] = load_te(path)
        else:
            log.warning("no TE file for %s (%s); piRNA support uses reproducibility only", org, path)
    return out


def unannotated_calls(table_val: pd.DataFrame, pred_cls) -> np.ndarray:
    pred = np.asarray(pred_cls, dtype=object)
    unann = np.array([l is None or pd.isna(l) for l in table_val["label"]], dtype=bool)
    carved = table_val["carved"].to_numpy(dtype=bool) if "carved" in table_val \
        else np.zeros(len(table_val), bool)
    called = np.array([p is not None and not pd.isna(p) and p != AMBIGUOUS for p in pred],
                      dtype=bool)
    return unann & ~carved & called


def supported(table_val: pd.DataFrame, pred_cls, te_index) -> np.ndarray:
    """Boolean per row: unannotated call that is supported. `te_index` is None, a TEIndex
    (applied to every row) or {organism: TEIndex} (organisms absent from it are judged by
    reproducibility alone)."""
    pred = np.asarray(pred_cls, dtype=object)
    calls = unannotated_calls(table_val, pred)
    rf = table_val["repro_frac"].to_numpy(float) if "repro_frac" in table_val \
        else np.full(len(table_val), np.nan)
    out = calls & (rf >= REPRO_MIN)          # NaN compares False
    todo = np.flatnonzero(out & (pred == "piRNA"))
    if te_index is None or not len(todo) or (isinstance(te_index, dict) and not te_index):
        return out
    orgs = table_val["organism"].to_numpy(dtype=object) if "organism" in table_val else None
    chrom, start, end = (table_val[c].to_numpy() for c in ("chrom", "start", "end"))
    for i in todo:
        te = te_index if isinstance(te_index, TEIndex) else \
            (te_index.get(orgs[i]) if orgs is not None else None)
        if te is not None and not te.overlaps(chrom[i], start[i], end[i]):
            out[i] = False
    return out


def supported_rate(table_val, pred_cls, te_index) -> tuple[float | None, int]:
    """(share of unannotated calls that are supported, number of unannotated calls)."""
    n = int(unannotated_calls(table_val, pred_cls).sum())
    if not n:
        return None, 0
    return float(supported(table_val, pred_cls, te_index).sum() / n), n
