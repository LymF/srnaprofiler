"""The biology gate as the benchmark sees it: gate_<class> columns of a clusters.tsv.

`classify` writes one 0/1 column per enabled class, computed by RuleScorer.passes_gate on
the record. The size evidence reads record.size_counts (per-length read mass), which the
TSV cannot carry: rebuilding it from the szp_*/szm_* fractions is not exact (both are
rounded to 4 decimals, and the strength rounds k = p_hat * n), so a rebuilt record can
gate differently from the original. The benchmark therefore reads the columns instead of
rebuilding records; production decides on the records themselves.
"""
import numpy as np

GATE_PREFIX = "gate_"


def gate_columns(table) -> list[str]:
    return [c for c in table.columns if str(c).startswith(GATE_PREFIX)]


def allowed_matrix(table, classes) -> np.ndarray:
    """(rows, classes) bool: does the row pass the gate of each class? A class with no
    column (not enabled in the profile) or a NaN cell passes nothing."""
    cols = gate_columns(table)
    if not cols:
        raise ValueError("the table has no gate_<class> columns: rebuild the benchmark "
                         "clusters.tsv with the current classify")
    out = np.zeros((len(table), len(classes)), bool)
    for j, c in enumerate(classes):
        col = f"{GATE_PREFIX}{c}"
        if col in table:
            v = table[col].to_numpy(float)
            out[:, j] = np.nan_to_num(v, nan=0.0) >= 0.5
    return out
