"""Ablation: rules + biology gate, no ML.

For every CV fold, the rules' class is kept when its gate_<class> column passes and becomes
ambiguous otherwise (miRNA and ambiguous rows are untouched). Written in the folds.jsonl
format as algorithm "rules_gate" so it can be compared with the ML rows of the same folds.

Run from the repo root:  python benchmark/scripts/ablate_gate.py <profile> <results_dir> <out_dir>
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, ".")
from benchmark.ml.cv import _te_indexes, folds_for, tuning_accessions
from benchmark.ml.dataset import build_table
from benchmark.ml.evaluate import AMBIGUOUS, fold_metrics
from mlgate import GATE_PREFIX


def gated_rules(table):
    rule = table["rule_cls"].to_numpy(dtype=object)
    out = rule.copy()
    for cls in set(rule) - {AMBIGUOUS, "miRNA"}:
        col = GATE_PREFIX + str(cls)
        rows = rule == cls
        if col not in table:
            out[rows] = AMBIGUOUS
            continue
        ok = np.nan_to_num(table[col].to_numpy(float), nan=0.0) >= 0.5
        out[rows & ~ok] = AMBIGUOUS
    return out


def main(profile, results_dir, out_dir):
    manifest = Path("benchmark/manifest.tsv")
    table = build_table(tuning_accessions(manifest, profile), Path(results_dir), manifest)
    table = table[table["profile"] == profile].reset_index(drop=True)
    te = _te_indexes(tuple(sorted(set(table["organism"]))))
    pred_all = gated_rules(table)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "folds.jsonl", "w") as fh:
        for unit, tr, va in folds_for(table, profile):
            val = table.loc[va]
            pred = pred_all[table.index.get_indexer(va)]
            m = fold_metrics(val, pred, np.zeros(len(val)), pred, None, te)
            fh.write(json.dumps({"profile": profile, "variant": "main", "algorithm": "rules_gate",
                                 "seed": 0, "unit": unit, "seeds": 1, "t": None,
                                 "n_train": int(len(tr)), "n_val": int(len(val)),
                                 "metrics": m}) + "\n")
    changed = float(np.mean(pred_all != table["rule_cls"].to_numpy(dtype=object)))
    print(f"{profile}: share of loci the gate turned ambiguous: {changed:.4f}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
