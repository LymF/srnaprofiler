"""Cross-validation driver: (variant, algorithm, seed, fold) -> one JSON line in folds.jsonl.

Resumable: lines already present (without "error") are skipped; each result is appended and
flushed by the parent process as soon as it finishes. A fold that raises writes a line with
"error" and no metrics (retried on the next run; ml-select refuses while it is unresolved).
Real models take fit(X, y) while evaluate.choose_threshold calls fit(table); `evaluate.TableFit`
adapts (y = table["label"]).
"""
import functools
import logging
import json
import shutil
import tempfile
from pathlib import Path

import joblib
import pandas as pd

from benchmark.manifest import load_manifest
from benchmark.ml.bioaudit import TRUTH_DIR, load_te_indexes
from benchmark.ml.dataset import build_table, load_ml_config
from benchmark.ml.evaluate import (TableFit, choose_threshold, complete_rules, decide,
                                   fold_metrics)
from benchmark.ml.folds import split
from mlmodels import REGISTRY
from mlfeatures import FEATURES, REPRO_COLS, RULE_COLS

log = logging.getLogger(__name__)
# "bio": only the textbook small-RNA signals (read size, strand, 1U/10A/5A, ping-pong, phasing);
# no abundance, mapping or locus-shape features, which can carry species-specific shortcuts.
_NOT_BIO = {"reads", "reads_rpm", "distinct_seqs", "multimap_frac", "unique_frac",
            "density_pvalue", "read_span", "start_entropy", "top_start_frac", "starts_per_kb"}
VARIANTS = {"main": list(FEATURES), "stack": list(FEATURES) + list(RULE_COLS),
            "repro": list(FEATURES) + list(REPRO_COLS),
            "bio": [f for f in FEATURES if f not in _NOT_BIO and not f.startswith("ctx_")]}
XSPECIES = ("Arabidopsis thaliana", "Oryza sativa")


def folds_for(table: pd.DataFrame, profile: str):
    """[(unit, train_index_labels, val_index_labels)] incl. plant cross-species folds."""
    out = [(unit, tr, va) for tr, va, unit in split(table)]
    if profile == "plant":
        a, b = XSPECIES
        ia, ib = table.index[table["organism"] == a], table.index[table["organism"] == b]
        if len(ia) and len(ib):
            out.append((f"xspecies:{b}", ia.to_numpy(), ib.to_numpy()))
            out.append((f"xspecies:{a}", ib.to_numpy(), ia.to_numpy()))
    return out


@functools.lru_cache(maxsize=4)
def _te_indexes(organisms: tuple) -> dict:
    return load_te_indexes(load_ml_config(), organisms=set(organisms))


@functools.lru_cache(maxsize=2)
def _load(path: str) -> pd.DataFrame:
    return joblib.load(path)


def _run_fold(table_path, profile, variant, algorithm, seed, unit, tr, va, nseeds=1):
    key = {"profile": profile, "variant": variant, "algorithm": algorithm, "seed": seed,
           "unit": unit, "seeds": nseeds}
    try:
        table = _load(table_path)
        train, val = table.loc[tr], table.loc[va]
        feats = VARIANTS[variant]
        rule = val["rule_cls"].to_numpy(dtype=object)
        if algorithm == "rules":
            model = REGISTRY["rules"](seed, feats)
            t, pred, conf = None, rule, None
        else:
            def factory():
                return TableFit(REGISTRY[algorithm](seed, feats))
            t = choose_threshold(factory, train, feats, seed)
            model = factory().fit(train)
            pred = decide(model, val, feats, t).to_numpy(dtype=object)
            proba = model.predict_proba(val)
            conf = proba.to_numpy().max(axis=1) if proba.shape[1] else None
        p_real = model.p_real(val)
        te = _te_indexes(tuple(sorted(set(table['organism']))))
        metrics = fold_metrics(val, pred, p_real, rule, conf, te)
        comp = complete_rules(rule, pred).to_numpy(dtype=object)
        cmetrics = fold_metrics(val, comp, p_real, rule, None, te)
        return {**key, "t": t, "n_train": int(len(train)), "n_val": int(len(val)),
                "metrics": metrics, "complete_rules_metrics": cmetrics}
    except Exception as exc:  # recorded, never silently dropped
        return {**key, "error": repr(exc)}


def _drop_partial_tail(path: Path) -> None:
    """Truncate to the end of the last complete line (a killed run may leave a partial one)."""
    if not path.exists():
        return
    data = path.read_bytes()
    if data and not data.endswith(b"\n"):
        with open(path, "r+b") as fh:
            fh.truncate(data.rfind(b"\n") + 1)


def read_lines(path: Path) -> list[dict]:
    path = Path(path)
    if not path.exists():
        return []
    out = []
    rows = [l for l in path.read_text().splitlines() if l.strip()]
    for i, l in enumerate(rows):
        try:
            out.append(json.loads(l))
        except json.JSONDecodeError:
            if i != len(rows) - 1:
                raise
            log.warning("%s: skipping truncated final line", path)
    return out


def _key(r):
    return (r["profile"], r["variant"], r["algorithm"], int(r["seed"]), r["unit"])


def run_cv(table: pd.DataFrame, profile: str, algorithms: list[str], variants: list[str],
           seeds: int, out: Path, n_jobs: int = 8) -> int:
    """Run all pending combinations; returns the number of lines appended."""
    if "repro_frac" not in table:
        raise ValueError("table lacks repro_frac (needed by the biology metric; build_table "
                         "computes it with repro_features)")
    _te_indexes(tuple(sorted(set(table["organism"]))))   # fail early on a missing TE file
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    path = out / "folds.jsonl"
    table = table.reset_index(drop=True)
    algs = [a for a in algorithms if a != "rules"]
    folds = folds_for(table, profile)
    combos = [("main", "rules", 0)]   # baseline: once per fold, seed-independent
    combos +=[(v, a, s) for v in variants for a in algs for s in range(seeds)]
    done = {_key(r) for r in read_lines(path) if "error" not in r}
    todo = [(v, a, s, u, tr, va) for v, a, s in combos for u, tr, va in folds
            if (profile, v, a, s, u) not in done]
    if not todo:
        return 0
    tmp = Path(tempfile.mkdtemp(prefix="cvtable"))
    n = 0
    try:
        tpath = str(tmp / "table.joblib")
        joblib.dump(table, tpath)
        gen = joblib.Parallel(n_jobs=n_jobs, return_as="generator_unordered")(
            joblib.delayed(_run_fold)(tpath, profile, v, a, s, u, tr, va, seeds)
            for v, a, s, u, tr, va in todo)
        with open(path, "a") as fh:
            _drop_partial_tail(path)
            for rec in gen:
                fh.write(json.dumps(rec) + "\n")
                fh.flush()
                n += 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return n


def tuning_accessions(manifest: Path, profile: str) -> list[str]:
    return [l.accession for l in load_manifest(Path(manifest))
            if l.role == "tuning" and l.taxon_group == profile]


def cv_main(profile, algorithms, variants, seeds, results_dir, out, manifest, n_jobs) -> int:
    names = sorted(REGISTRY) if algorithms == ["all"] else algorithms
    bad = [a for a in names if a not in REGISTRY] + [v for v in variants if v not in VARIANTS]
    if bad:
        raise ValueError(f"unknown algorithm/variant: {bad}")
    table = build_table(tuning_accessions(manifest, profile), results_dir, manifest)
    table = table[table["profile"] == profile]
    return run_cv(table, profile, names, variants, seeds, out, n_jobs)
