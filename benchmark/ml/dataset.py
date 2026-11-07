"""Locus feature table from benchmark runs: one row per locus, labelled by truth match."""
from pathlib import Path

import pandas as pd

from benchmark.manifest import load_manifest
from benchmark.metrics_class import match_loci, read_loci
from benchmark.orchestrate import load_truth
from mlfeatures import FEATURES, REPRO_COLS, RULE_COLS, derive

__all__ = ["FEATURES", "RULE_COLS", "REPRO_COLS", "build_table"]

MIN_OVERLAP_FRAC = 0.5
META = ["acc", "organism", "tissue", "replicate", "profile", "chrom", "start", "end",
        "rule_cls", "label", "is_control", "carved"]


def _opt(value):
    return float("nan") if value in (None, "", "NA") else float(value)


def build_table(accessions: list[str], results_dir: Path, manifest: Path,
                final_test: bool = False) -> pd.DataFrame:
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
        pairs = match_loci(loci, load_truth(lib), MIN_OVERLAP_FRAC)
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
    for c in ("label", "rule_cls"):
        df[c] = df[c].astype(object).where(df[c].notna(), None)
    return df
