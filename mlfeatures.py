"""Locus feature derivation shared by the benchmark and the ML classifier mode.

Single source of feature computation: maps one locus (a clusters.tsv row or the
attributes of a ClusterRecord) to the numeric FEATURES vector. numpy/math only.
"""
import math

_NUMERIC = [
    "reads", "reads_rpm", "distinct_seqs", "size_purity", "frac_in_window", "frac_top",
    "pct_1u", "pct_10a", "pct_5a", "gc", "pingpong_z", "pingpong_pairs10", "pingpong_frac",
    "overlap_dup_frac", "phase_score", "phase_pvalue", "phase_number", "phase_ratio",
    "multimap_frac", "unique_frac", "density_pvalue", "read_span",
]
_SHAPE = ["start_entropy", "top_start_frac", "starts_per_kb"]
_CTX = ["ctx_te", "ctx_trna", "ctx_rrna", "ctx_gene"]

FEATURES = [*_NUMERIC, "dicer_len", *[f"szp_{L}" for L in range(18, 36)],
            *[f"szm_{L}" for L in range(18, 36)], *_SHAPE, *_CTX]
# Note: with `classify --annotation` the rule `confidence` is nudged by the annotation
# (and `review` set), so it must never be used to train or evaluate an ML model.
# Benchmarks and training runs are made without --annotation.
RULE_COLS = ["score", "confidence", "s_1u", "s_5a", "s_fractop", "s_fractop_bal",
             "s_fractop_bid", "s_pingpong", "s_phase"]
REPRO_COLS = ["repro_frac", "repro_rpm_ratio"]

_NAN = float("nan")


def _num(value) -> float:
    if value is None:
        return _NAN
    if isinstance(value, str):
        v = value.strip()
        if v in ("", "NA", "N", "nan", "NaN"):
            return _NAN
        try:
            return float(v)
        except ValueError:
            return _NAN
    try:
        return float(value)
    except (TypeError, ValueError):
        return _NAN


def derive(row) -> dict[str, float]:
    """FEATURES values for one locus; "NA"/None/"" become NaN, ctx_* are NaN for now."""
    out = {}
    for name in FEATURES:
        if name in _CTX:
            out[name] = _NAN
        elif name == "dicer_len":
            out[name] = _num(row.get("dicer_call"))
        else:
            out[name] = _num(row.get(name))
    return out
