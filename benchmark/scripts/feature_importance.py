"""Feature-group importance of one ML algorithm on the CV folds (tuning data only).

Per leave-one-unit-out fold: fit the model, then (a) permute each feature group on the
validation rows (rows shuffled jointly inside the group) and record the metric change with
the model fixed, and (b) for the groups in RETRAIN refit without the group. The biology gate
columns are not features, so they stay intact: the change isolates the model.

Run from the repo root:
  python benchmark/scripts/feature_importance.py <profile> <algorithm> <results_dir> <out.jsonl> [n_jobs]
"""
import json
import sys
from pathlib import Path

import joblib
import numpy as np

sys.path.insert(0, ".")
from benchmark.ml.cv import VARIANTS, _te_indexes, folds_for, tuning_accessions
from benchmark.ml.dataset import build_table
from benchmark.ml.evaluate import TableFit, decide, fold_metrics
from mlmodels import REGISTRY

GROUPS = {
    "size": lambda f: f.startswith(("szp_", "szm_")) or f in ("size_purity", "dicer_len", "frac_in_window"),
    "strand": lambda f: f == "frac_top",
    "nucleotide": lambda f: f in ("pct_1u", "pct_10a", "pct_5a", "gc"),
    "pingpong": lambda f: f.startswith("pingpong") or f == "overlap_dup_frac",
    "phasing": lambda f: f.startswith("phase_"),
    "abundance": lambda f: f in ("reads", "reads_rpm", "distinct_seqs"),
    "mapping": lambda f: f in ("multimap_frac", "unique_frac"),
    "shape": lambda f: f in ("density_pvalue", "read_span", "start_entropy", "top_start_frac", "starts_per_kb"),
    "context": lambda f: f.startswith("ctx_"),
}
RETRAIN = ("context", "abundance")
KEYS = ("macro_f1", "control_real_rate", "unann_supported_rate", "ambiguous_share")
REPEATS = 3


def score(model, val, feats, te):
    pred = decide(model, val, feats, 0.0).to_numpy(dtype=object)
    m = fold_metrics(val, pred, model.p_real(val), val["rule_cls"].to_numpy(dtype=object), None, te)
    return {k: m[k] for k in KEYS}


def one_fold(table, profile, algorithm, unit, tr, va, te):
    feats = list(VARIANTS["main"])
    train, val = table.loc[tr], table.loc[va]
    model = TableFit(REGISTRY[algorithm](0, feats)).fit(train)
    out = [{"unit": unit, "group": "none", "kind": "base", **score(model, val, feats, te)}]
    rng = np.random.default_rng(0)
    for g, member in GROUPS.items():
        cols = [f for f in feats if member(f)]
        if not cols:
            continue
        for rep in range(REPEATS):
            shuffled = val.copy()
            perm = rng.permutation(len(val))
            shuffled[cols] = val[cols].to_numpy()[perm]
            out.append({"unit": unit, "group": g, "kind": "permute", "rep": rep,
                        **score(model, shuffled, feats, te)})
        if g in RETRAIN:
            rest = [f for f in feats if f not in cols]
            m2 = TableFit(REGISTRY[algorithm](0, rest)).fit(train)
            out.append({"unit": unit, "group": g, "kind": "retrain",
                        **score(m2, val, rest, te)})
    return out


def main(profile, algorithm, results_dir, out, n_jobs=6):
    manifest = Path("benchmark/manifest.tsv")
    table = build_table(tuning_accessions(manifest, profile), Path(results_dir), manifest)
    table = table[table["profile"] == profile].reset_index(drop=True)
    te = _te_indexes(tuple(sorted(set(table["organism"]))))
    folds = [f for f in folds_for(table, profile) if not str(f[0]).startswith("xspecies")]
    res = joblib.Parallel(n_jobs=int(n_jobs))(
        joblib.delayed(one_fold)(table, profile, algorithm, u, tr, va, te) for u, tr, va in folds)
    with open(out, "w") as fh:
        for rows in res:
            for r in rows:
                fh.write(json.dumps(r) + "\n")
    print(f"{out}: {sum(len(r) for r in res)} lines")


if __name__ == "__main__":
    main(*sys.argv[1:6])
