#!/usr/bin/env python3
"""Regenerate the validation tables from the runs listed in a manifest.

Writes one tidy CSV per analysis, ready to plot. Numbers quoted in the roadmap
or in a manuscript should come from here, not from prose: rerun this and the
tables rebuild from the .clusters.tsv files themselves.

Usage:
    python3 extras/collect_validation.py docs/validation/manifest.tsv outdir/

Manifest columns (tab-separated, '#' comments allowed):
    run_id  kingdom  path  mincov  config  mirna_source  note
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_piclusters import load_reference, mark_inside, ARM_LEN
from validate_plant_tas import load_tas_loci

# TAIR10.1 chromosome lengths (RefSeq accessions)
PLANT_LEN = {"NC_003070.9": 30427671, "NC_003071.7": 19698289,
             "NC_003074.8": 23459830, "NC_003075.7": 18585056,
             "NC_003076.8": 26975502}

HERE = os.path.dirname(os.path.abspath(__file__))
PIC_BED = os.path.join(HERE, "piC_ref_dm6.bed")


def _load(path):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d["is_review"] = d.review.astype(str).str.lower() == "true"
    return d


def summary_table(runs):
    """One row per run: cluster counts, class shares, review rate."""
    rows = []
    for r in runs.itertuples(index=False):
        d = _load(r.path)
        counts = d.cls.value_counts()
        row = dict(run_id=r.run_id, kingdom=r.kingdom, mincov=r.mincov,
                   config=r.config, mirna_source=r.mirna_source,
                   clusters=len(d), review_rate=round(d.is_review.mean(), 4),
                   kept=int((~d.is_review).sum()),
                   reads_total=round(d.reads.sum(), 1))
        for cls, n in counts.items():
            row[f"n_{cls}"] = int(n)
        rows.append(row)
    return pd.DataFrame(rows)


def size_distribution(runs):
    """DicerCall shares per run — the size profile that separates the kingdoms."""
    rows = []
    for r in runs.itertuples(index=False):
        d = _load(r.path)
        dc = pd.to_numeric(d.dicer_call, errors="coerce")
        for size, frac in dc.value_counts(normalize=True).items():
            if not np.isnan(size):
                rows.append(dict(run_id=r.run_id, kingdom=r.kingdom,
                                 size=int(size), fraction=round(frac, 5)))
    return pd.DataFrame(rows)


def pic_recall(runs):
    """Drosophila: recall over known piRNA clusters and read-mass enrichment."""
    ref = load_reference(PIC_BED)
    genome = sum(ARM_LEN.values())
    ref_bp = sum(min(e, ARM_LEN[c]) - s for c, s, e in ref[["chrom", "start", "end"]].values)
    rows = []
    for r in runs[runs.kingdom == "arthropod"].itertuples(index=False):
        d = _load(r.path)
        d = d[d.chrom.isin(ARM_LEN)]
        if d.empty:
            continue
        d = d.assign(inside=mark_inside(d, ref))
        pir = d[(d.cls == "piRNA") & ~d.is_review]
        rec = 0
        for chrom, g in ref.groupby("chrom"):
            p = pir[pir.chrom == chrom]
            for rs, re_ in g[["start", "end"]].values:
                if ((p.start < re_) & (p.end > rs)).any():
                    rec += 1
        row = dict(run_id=r.run_id, n_reference=len(ref), recovered=rec,
                   recall=round(rec / len(ref), 4))
        for cls in ("piRNA", "siRNA"):
            m = d.cls == cls
            if m.any():
                share = d.reads[m & d.inside].sum() / d.reads[m].sum()
                row[f"{cls}_enrichment"] = round(share / (ref_bp / genome), 3)
        rows.append(row)
    return pd.DataFrame(rows)


def tas_recall(runs, gff):
    """Plants: are the TAS loci called phasiRNA?"""
    tas = load_tas_loci(gff)
    rows = []
    for r in runs[runs.kingdom == "plant"].itertuples(index=False):
        d = _load(r.path)
        for t in tas.itertuples(index=False):
            g = d[(d.chrom == t.chrom) & (d.start < t.end) & (d.end > t.start)]
            if g.empty:
                rows.append(dict(run_id=r.run_id, locus=t.name, called=False))
                continue
            top = g.loc[g.reads.idxmax()]
            rows.append(dict(run_id=r.run_id, locus=t.name, called=True,
                             cls=top.cls, dicer_call=top.dicer_call,
                             reads=round(top.reads, 1),
                             phase_score=top.phase_score,
                             phase_pvalue=top.phase_pvalue))
    return pd.DataFrame(rows)


def te_enrichment(runs, te_path):
    """Plants: class vs transposon overlap (the RdDM prediction)."""
    te = pd.read_csv(te_path, sep="\t")
    chrom_of = {"1": "NC_003070.9", "2": "NC_003071.7", "3": "NC_003074.8",
                "4": "NC_003075.7", "5": "NC_003076.8"}
    te["chrom"] = te.Transposon_Name.str.extract(r"AT(\d)TE")[0].map(chrom_of)
    te = te.dropna(subset=["chrom"])
    masks = {}
    for c, g in te.groupby("chrom"):
        m = np.zeros(PLANT_LEN[c] + 1, bool)
        for s, e in g[["Transposon_min_Start", "Transposon_max_End"]].values:
            m[int(s):int(e)] = True
        masks[c] = m
    baseline = sum(m.sum() for m in masks.values()) / sum(PLANT_LEN.values())
    rows = []
    for r in runs[runs.kingdom == "plant"].itertuples(index=False):
        d = _load(r.path)
        d = d[d.chrom.isin(PLANT_LEN) & ~d.is_review]
        if d.empty:
            continue
        inside = np.zeros(len(d), bool)
        for c, m in masks.items():
            for k in np.where((d.chrom == c).values)[0]:
                s, e = d.start.values[k], min(d.end.values[k], PLANT_LEN[c])
                if e > s:
                    inside[k] = m[s:e].mean() > 0.5
        d = d.assign(in_te=inside)
        for cls, g in d.groupby("cls"):
            rows.append(dict(run_id=r.run_id, cls=cls, n=len(g),
                             frac_in_te=round(g.in_te.mean(), 4),
                             enrichment=round(g.in_te.mean() / baseline, 3),
                             genome_baseline=round(baseline, 4)))
    return pd.DataFrame(rows)


def main(manifest_path, outdir):
    runs = pd.read_csv(manifest_path, sep="\t", comment="#")
    runs = runs[runs.path.apply(os.path.exists)]
    if runs.empty:
        sys.exit("No run in the manifest points at an existing file.")
    os.makedirs(outdir, exist_ok=True)
    cfg = dict(gff=os.environ.get("TAIR_GFF"), te=os.environ.get("TAIR_TE"))

    tables = {"summary": summary_table(runs),
              "size_distribution": size_distribution(runs),
              "pic_recall": pic_recall(runs)}
    if cfg["gff"] and os.path.exists(cfg["gff"]):
        tables["tas_recall"] = tas_recall(runs, cfg["gff"])
    if cfg["te"] and os.path.exists(cfg["te"]):
        tables["te_enrichment"] = te_enrichment(runs, cfg["te"])

    for name, t in tables.items():
        if t is None or t.empty:
            continue
        path = os.path.join(outdir, f"{name}.csv")
        t.to_csv(path, index=False)
        print(f"{path}: {len(t)} rows")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
