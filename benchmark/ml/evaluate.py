"""Abstention threshold and per-fold evaluation of the ML classifiers.

fold_metrics scores one validation fold at locus level: the table already carries the truth
label of every locus, so precision/recall are computed directly from (pred, label) rows.
Precision matches the benchmark's classification_metrics (abstentions are excluded from tp/fp,
unlabeled loci are never counted as false). Recall differs: here it is locus-level (labelled
loci of the class predicted as the class / labelled loci of the class), whereas the benchmark
reports region-level recall from real runs.
"""
import numpy as np
import pandas as pd

from benchmark.metrics_class import _ece
from benchmark.ml.folds import chromosome_units

AMBIGUOUS = "ambiguous"
THRESHOLDS = [round(0.05 * i, 2) for i in range(20)]
FALLBACK_T = 0.95
INNER_FOLDS = 3


def _macro_precision(pred, label):
    """Mean precision over classes with >= 1 annotated row predicted as that class
    (abstentions and unlabeled rows ignored). None when no class qualifies."""
    pred = np.asarray(pred, dtype=object)
    label = np.asarray(label, dtype=object)
    keep = np.array([p != AMBIGUOUS and p is not None and l is not None
                     for p, l in zip(pred, label)], dtype=bool)
    if not keep.any():
        return None
    p, l = pred[keep], label[keep]
    precs = [float(np.mean(l[p == c] == c)) for c in np.unique(p)]
    return float(np.mean(precs))


def _inner_groups(table, seed):
    """Row-position arrays, one per inner validation group (chromosome units grouped,
    seeded). With a single unit: seeded random row split into 3 (documented fallback)."""
    rng = np.random.default_rng(seed)
    units = chromosome_units(table).to_numpy()
    names = np.array(sorted(set(units)))
    if len(names) < 2:
        order = rng.permutation(len(table))
        return [np.sort(order[i::INNER_FOLDS]) for i in range(INNER_FOLDS)]
    k = min(INNER_FOLDS, len(names))
    names = rng.permutation(names)
    return [np.flatnonzero(np.isin(units, names[i::k])) for i in range(k)]


def choose_threshold(factory, table_train: pd.DataFrame, features: list[str], seed: int) -> float:
    """Smallest t in {0.00, 0.05, ..., 0.95} whose macro precision on annotated loci (max
    probability >= t) is >= the rules' macro precision on the same loci; 0.95 if none.

    `factory` is a zero-arg callable returning a fresh unfitted Model. Probabilities are pooled
    from an inner 3-fold split by chromosome unit (units grouped, seeded; as many groups as
    units, minimum 2). With a single unit the split falls back to a seeded random 3-fold split
    of rows. Carved rows are not scored. The rules are compared on the same annotated rows.
    """
    table_train = table_train.reset_index(drop=True)
    groups = _inner_groups(table_train, seed)
    conf, pred, idx = [], [], []
    for g in groups:
        val_mask = np.zeros(len(table_train), bool)
        val_mask[g] = True
        if val_mask.all() or not val_mask.any():
            continue
        model = factory().fit(table_train[~val_mask])
        val = table_train[val_mask]
        proba = model.predict_proba(val)
        arr = proba.to_numpy()
        conf.append(arr.max(axis=1))
        pred.append(np.asarray(proba.columns, dtype=object)[arr.argmax(axis=1)])
        idx.append(np.flatnonzero(val_mask))
    if not conf:
        return FALLBACK_T
    conf, pred, idx = np.concatenate(conf), np.concatenate(pred), np.concatenate(idx)
    label = table_train["label"].to_numpy(dtype=object)[idx]
    rule = table_train["rule_cls"].to_numpy(dtype=object)[idx]
    carved = table_train["carved"].to_numpy(dtype=bool)[idx] if "carved" in table_train \
        else np.zeros(len(idx), bool)
    ann = np.array([l is not None and not pd.isna(l) for l in label]) & ~carved
    conf, pred, label, rule = conf[ann], pred[ann], label[ann], rule[ann]
    target = _macro_precision(rule, label)
    target = 0.0 if target is None else target
    for t in THRESHOLDS:
        shown = np.where(conf >= t, pred, AMBIGUOUS)
        mp = _macro_precision(shown, label)
        if mp is not None and mp >= target:
            return t
    return FALLBACK_T


def decide(model, table: pd.DataFrame, features: list[str], t: float) -> pd.Series:
    """Argmax class, or "ambiguous" when the max probability is < t; carved rows keep miRNA."""
    proba = model.predict_proba(table)
    arr = proba.to_numpy()
    cls = np.asarray(proba.columns, dtype=object)[arr.argmax(axis=1)]
    out = np.where(arr.max(axis=1) < t, AMBIGUOUS, cls).astype(object)
    if "carved" in table:
        out[table["carved"].to_numpy(dtype=bool)] = "miRNA"
    return pd.Series(out, index=table.index, dtype=object)


def complete_rules(rule_cls, pred_cls) -> pd.Series:
    """Keep the rule class where it is not ambiguous, take the ML class elsewhere."""
    rule = pd.Series(rule_cls)
    pred = pd.Series(pred_cls, index=rule.index)
    return rule.where(rule != AMBIGUOUS, pred)


def _per_class(pred, label):
    classes = sorted({c for c in list(pred) + list(label)
                      if c is not None and not pd.isna(c) and c != AMBIGUOUS})
    out = {}
    for c in classes:
        tp = sum(p == c and l == c for p, l in zip(pred, label))
        fp = sum(p == c and l is not None and l != c for p, l in zip(pred, label))
        n_true = sum(l == c for l in label)
        out[c] = {"tp": tp, "fp": fp, "n_truth": n_true,
                  "precision": tp / (tp + fp) if tp + fp else None,
                  "recall": tp / n_true if n_true else None}
    return out


def fold_metrics(table_val: pd.DataFrame, pred_cls, p_real, rule_cls, confidence=None) -> dict:
    """Metrics of one validation fold. `confidence` (max class probability per row) feeds the
    ECE over annotated, non-abstained rows; ece is None without it. See module docstring for
    the locus-level recall caveat."""
    pred = np.asarray(pred_cls, dtype=object)
    label = np.array([None if (l is None or pd.isna(l)) else l for l in table_val["label"]],
                     dtype=object)
    p_real = np.asarray(p_real, float)
    rule = np.asarray(rule_cls, dtype=object)
    per_class = _per_class(pred, label)
    f1s = []
    for c, m in per_class.items():
        if m["n_truth"] == 0:
            continue
        p, r = m["precision"] or 0.0, m["recall"] or 0.0
        f1s.append(2 * p * r / (p + r) if p + r else 0.0)
    gap = None
    if "repro_frac" in table_val:
        rf = table_val["repro_frac"].to_numpy(float)
        hi, lo = rf[(p_real >= 0.5) & ~np.isnan(rf)], rf[(p_real < 0.5) & ~np.isnan(rf)]
        if len(hi) and len(lo):
            gap = float(hi.mean() - lo.mean())
    ctrl = table_val["is_control"].to_numpy(dtype=bool) if "is_control" in table_val \
        else np.zeros(len(table_val), bool)
    ece = None
    if confidence is not None:
        conf = np.asarray(confidence, float)
        ok = np.array([l is not None and p != AMBIGUOUS for p, l in zip(pred, label)], bool)
        ece = float(_ece([(float(c), bool(p == l)) for c, p, l
                          in zip(conf[ok], pred[ok], label[ok])]))
    return {
        "per_class": per_class,
        "macro_f1": float(np.mean(f1s)) if f1s else None,
        "ambiguous_share": float(np.mean(pred == AMBIGUOUS)) if len(pred) else None,
        "real_repro_gap": gap,
        "control_real_rate": float(np.mean(p_real[ctrl] >= 0.5)) if ctrl.any() else None,
        "rules_control_rate": float(np.mean(rule[ctrl] != AMBIGUOUS)) if ctrl.any() else None,
        "ece": ece,
    }
