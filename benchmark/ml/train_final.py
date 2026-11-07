"""Freeze a model: fit on every tuning row of a profile, choose t, write joblib + schema.json."""
import json
import subprocess
from pathlib import Path

import joblib
import pandas as pd

from benchmark.manifest import load_manifest
from benchmark.ml.cv import VARIANTS, tuning_accessions
from benchmark.ml.dataset import build_table
from benchmark.ml.evaluate import TableFit, choose_threshold
from mlmodels import REGISTRY
from benchmark.run import REPO, pipeline_sha256
from mlfeatures import FEATURES


def _commit() -> str:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                             text=True, check=True)
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def train_and_save(table: pd.DataFrame, profile: str, algorithm: str, variant: str,
                   seed: int, out: Path) -> dict:
    if algorithm not in REGISTRY or algorithm == "rules":
        raise ValueError(f"unknown algorithm {algorithm!r}")
    if variant != "main":
        raise ValueError("only the main variant can be shipped (it must not use rule outputs)")
    feats = VARIANTS[variant]
    if feats != list(FEATURES):
        raise ValueError("main variant features differ from mlfeatures.FEATURES")
    table = table.reset_index(drop=True)
    t = choose_threshold(lambda: TableFit(REGISTRY[algorithm](seed, feats)), table, feats, seed)
    model = REGISTRY[algorithm](seed, feats).fit(table, table["label"])
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(model, out)
    schema = {"features": list(FEATURES), "classes": [str(c) for c in model.classes_],
              "threshold": float(t), "algorithm": algorithm, "variant": variant,
              "seed": int(seed), "train_commit": _commit(), "pipeline_hash": pipeline_sha256()}
    out.with_suffix(".schema.json").write_text(json.dumps(schema, indent=2) + "\n")
    return schema


def train_final(profile, algorithm, variant, seed, out, results_dir, manifest) -> dict:
    accs = tuning_accessions(Path(manifest), profile)
    if not accs:
        raise ValueError(f"no tuning library for profile {profile!r} in {manifest}")
    table = build_table(accs, results_dir, manifest)
    table = table[table["profile"] == profile]
    return train_and_save(table, profile, algorithm, variant, seed, out)
