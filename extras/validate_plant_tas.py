#!/usr/bin/env python3
"""Validate a plant run against the canonical Arabidopsis TAS loci.

The TAS genes (TAS1a/b/c, TAS2, TAS3a/b, TAS4) are the textbook tasiRNA sources:
21-nt, strongly phased, triggered by miR173/miR390. They are to a plant run what
flamenco and 42AB are to the Drosophila one -- a landmark whose expected
signature is published, and which nothing in the pipeline knows about.

Coordinates are read straight out of a TAIR10 GFF3 by gene name, so nothing is
hard-coded here and the reference cannot silently drift from the assembly the
reads were mapped to. Get the GFF3 from https://www.arabidopsis.org/download/
(Genes/TAIR10_genome_release) or Ensembl Plants.

Usage:
    python3 extras/validate_plant_tas.py run.clusters.tsv TAIR10_GFF3_genes.gff

Expected if the plant config works: the TAS loci come out phasiRNA, 21-nt,
with a significant phase_pvalue and high pct_1u (AGO1 loading); genome-wide,
24-nt hc-siRNA loci should dominate and carry high pct_5a (AGO4).
"""
import re
import sys

import numpy as np
import pandas as pd

TAS_NAME = re.compile(r"(?:^|[;=])(?:Name|Alias)=(TAS\d[a-z]?)\b", re.I)


def load_tas_loci(gff_path):
    """Gene features whose Name/Alias looks like TAS1a, TAS2, TAS3a..."""
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "transcribed_nonprotein_coding_gene"):
                continue
            m = TAS_NAME.search(f[8])
            if m:
                rows.append((m.group(1).upper(), f[0], int(f[3]) - 1, int(f[4])))
    df = pd.DataFrame(rows, columns=["name", "chrom", "start", "end"])
    return df.drop_duplicates(subset="name").sort_values("name")


def main(tsv_path, gff_path):
    tas = load_tas_loci(gff_path)
    if tas.empty:
        sys.exit("No TAS gene found in the GFF3 — check that it is a TAIR10 gene "
                 "annotation with Name= attributes.")
    d = pd.read_csv(tsv_path, sep="\t", low_memory=False)

    # The GFF3 uses 1/2/3..., a BAM may use Chr1/chr1: match on the trailing digits.
    def norm(c):
        return str(c).lower().replace("chr", "")
    d = d.assign(cnorm=d.chrom.map(norm))
    tas = tas.assign(cnorm=tas.chrom.map(norm))
    if not set(tas.cnorm) & set(d.cnorm):
        sys.exit(f"No chromosome in common. GFF3: {sorted(set(tas.cnorm))[:5]} "
                 f"vs TSV: {sorted(set(d.cnorm))[:5]}")

    print(f"{len(tas)} TAS loci found in the annotation\n")
    print(f"{'locus':8s} {'loci':>5s} {'reads':>10s} {'class':>10s} {'size':>5s} "
          f"{'phase_p':>9s} {'1U':>5s}")
    hits = 0
    for t in tas.itertuples(index=False):
        g = d[(d.cnorm == t.cnorm) & (d.start < t.end) & (d.end > t.start)]
        if g.empty:
            print(f"{t.name:8s} {'-':>5s} {'no locus called':>30s}")
            continue
        hits += 1
        top = g.loc[g.reads.idxmax()]
        w = g.reads / g.reads.sum()
        print(f"{t.name:8s} {len(g):5d} {g.reads.sum():10.0f} {top.cls:>10s} "
              f"{str(top.dicer_call):>5s} {g.phase_pvalue.min():9.2e} "
              f"{np.average(g.pct_1u, weights=w):5.3f}")
    print(f"\nTAS loci with at least one called locus: {hits}/{len(tas)}")

    dc = pd.to_numeric(d.dicer_call, errors="coerce")
    print("\ngenome-wide sanity (plants should be 24-nt dominated):")
    print("   size shares:", {int(k): round(v, 3)
                              for k, v in dc.value_counts(normalize=True).head(5).items()})
    print("   classes:", {k: round(v, 3) for k, v in d.cls.value_counts(normalize=True).items()})
    if "pct_5a" in d:
        for lo, hi, lbl in ((21, 22, "21-22nt"), (23, 24, "23-24nt")):
            m = (dc >= lo) & (dc <= hi)
            if m.any():
                print(f"   {lbl}: pct_1u={d.pct_1u[m].mean():.3f}  pct_5a={d.pct_5a[m].mean():.3f}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
