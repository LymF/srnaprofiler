"""Locus feature table from benchmark runs: one row per locus, labelled by truth match."""
import logging
from pathlib import Path

import pandas as pd
import yaml

from benchmark.manifest import load_manifest
from benchmark.metrics_class import match_loci, read_loci
from benchmark.orchestrate import load_truth
from benchmark.ml.controls import is_control
from benchmark.ml.replicates import repro_features
from mlfeatures import FEATURES, REPRO_COLS, RULE_COLS, derive

__all__ = ["FEATURES", "RULE_COLS", "REPRO_COLS", "build_table"]

CONFIG = Path(__file__).with_name("config.yaml")
TRUTH_DIR = Path(__file__).resolve().parents[1] / "truth"

log = logging.getLogger(__name__)


def load_ml_config(path: Path = CONFIG) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


META = ["acc", "organism", "tissue", "replicate", "profile", "chrom", "start", "end",
        "rule_cls", "label", "is_control", "carved"]


def _opt(value):
    return float("nan") if value in (None, "", "NA") else float(value)


def build_table(accessions: list[str], results_dir: Path, manifest: Path,
                final_test: bool = False) -> pd.DataFrame:
    cfg = load_ml_config()
    min_overlap_frac = float(cfg["min_overlap_frac"])
    libs = {l.accession: l for l in load_manifest(Path(manifest))}
    chosen = []
    for acc in accessions:
        if acc not in libs:
            raise ValueError(f"unknown accession {acc} in {manifest}")
        if libs[acc].role == "test" and not final_test:
            raise ValueError(f"{acc} is a test library; pass final_test=True to use it")
        chosen.append(libs[acc])
    rows = []
    for lib in chosen:
        path = Path(results_dir) / lib.accession / "work" / f"{lib.accession}.clusters.tsv"
        loci = read_loci(path)
        pairs = match_loci(loci, load_truth(lib), min_overlap_frac)
        for loc, region in pairs:
            rule_cls = loc.get("cls")
            row = {"acc": lib.accession, "organism": lib.organism, "tissue": lib.tissue,
                   "replicate": lib.replicate, "profile": lib.taxon_group,
                   "chrom": loc["chrom"], "start": loc["start"], "end": loc["end"],
                   "rule_cls": rule_cls, "label": region.cls if region else None,
                   "is_control": False, "carved": rule_cls == "miRNA"}
            row.update(derive(loc))
            for c in RULE_COLS:
                row[c] = _opt(loc.get(c))
            for c in REPRO_COLS:
                row[c] = float("nan")
            rows.append(row)
    df = pd.DataFrame(rows, columns=[*META, *FEATURES, *RULE_COLS, *REPRO_COLS])
    df = df.reset_index(drop=True)   # folds.split yields index labels: keep them unique
    df["is_control"] = False
    for org in df["organism"].unique():
        sub = cfg.get("controls", {}).get(org)
        bed = TRUTH_DIR / sub / "controls.bed" if sub else None
        if bed is None or not bed.exists():
            log.warning("no controls.bed for %s (%s); is_control stays False", org, bed)
            continue
        m = df["organism"] == org
        df.loc[m, "is_control"] = is_control(df[m], bed, min_overlap_frac)
    if len(df):
        df[REPRO_COLS] = repro_features(df, min_overlap_frac)[REPRO_COLS]
    for c in ("label", "rule_cls"):
        df[c] = df[c].astype(object).where(df[c].notna(), None)
    return df
