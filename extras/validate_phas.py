#!/usr/bin/env python3
"""Validate a plant run against annotated PHAS loci.

The six TAS genes are too narrow a reference to calibrate phasing on: two of
them are barely expressed in flower buds, so 4 loci end up carrying the whole
judgement. This checks against the PHAS annotation of sRNAanno instead, which
gives 35 loci for Arabidopsis and hundreds for the grasses, and reports:

  1. recall     -- PHAS loci overlapped by one of our phasiRNA calls
  2. phasing    -- our phase_pvalue inside the PHAS loci against the background,
                   which is the point: the test should separate them
  3. agreement  -- our p against the reference's own, locus by locus

Reference: Chen et al. (2021), Horticulture Research 8:45,
https://doi.org/10.1038/s41438-021-00480-8 (http://www.plantsrnas.org)

Usage:  python3 extras/validate_phas.py run.clusters.tsv [phas.gff3] [--cls phasiRNA]
"""
import re
import sys

import numpy as np
import pandas as pd

DEFAULT_GFF = "extras/phas21_ref_at_TAIR10.gff3"
ATTR = re.compile(r"(\w+)=([^;]+)")


def load_phas(gff_path):
    """PHAS features with the reference's own phasing statistics."""
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or not f[2].startswith("PHAS"):
                continue
            a = dict(ATTR.findall(f[8]))
            rows.append((a.get("ID", ""), f[0], int(f[3]) - 1, int(f[4]),
                         float(a.get("Pvalue", "nan")),
                         float(a.get("Phas_Score", "nan"))))
    return pd.DataFrame(rows, columns=["id", "chrom", "start", "end",
                                       "ref_pvalue", "ref_score"])


def norm_chrom(c):
    """Chr1 / chr1 / 1 / NC_003070.9 -> a comparable key.

    RefSeq accessions carry no chromosome number, so TAIR10's five are mapped
    explicitly; anything else falls back to its trailing digits.
    """
    refseq = {"NC_003070": "1", "NC_003071": "2", "NC_003074": "3",
              "NC_003075": "4", "NC_003076": "5",
              "NC_000932": "pt", "NC_037304": "mt"}
    s = str(c)
    base = s.split(".")[0]
    if base in refseq:
        return refseq[base]
    m = re.search(r"(\d+)$", s.lower().replace("chr", ""))
    return m.group(1) if m else s.lower()


def main(tsv_path, gff_path, cls_name):
    phas = load_phas(gff_path)
    if phas.empty:
        sys.exit(f"No PHAS feature in {gff_path}")
    d = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    for col in ("phase_pvalue", "cls", "review"):
        if col not in d.columns:
            sys.exit(f"{tsv_path} has no '{col}' column — is this a classify run?")

    d = d.assign(k=d.chrom.map(norm_chrom))
    phas = phas.assign(k=phas.chrom.map(norm_chrom))
    if not set(phas.k) & set(d.k):
        sys.exit(f"No chromosome in common. Reference: {sorted(set(phas.k))[:5]} "
                 f"Run: {sorted(set(d.k))[:5]}")

    review = d.review.astype(str).str.lower() == "true"
    called = d[(d.cls == cls_name) & ~review.values]
    print(f"{len(phas)} PHAS loci in the reference, {len(called):,} non-review "
          f"{cls_name} loci in the run\n")

    # 1. recall, and our p-value at each reference locus
    hits, ours = 0, []
    for r in phas.itertuples():
        same = d[d.k == r.k]
        ov = same[(same.start < r.end) & (same.end > r.start)]
        if len(ov):
            best = ov.loc[ov.phase_pvalue.idxmin()]
            ours.append((r.id, r.ref_pvalue, best.phase_pvalue, best.cls,
                         str(best.review), best.reads))
        if ((ov.cls == cls_name) & (ov.review.astype(str).str.lower() == "false")).any():
            hits += 1
    print(f"1. RECALL: {hits}/{len(phas)} PHAS loci hold >=1 non-review {cls_name} locus")

    # 2. does our phasing test separate PHAS regions from the rest?
    inside = np.zeros(len(d), bool)
    for r in phas.itertuples():
        inside |= ((d.k == r.k) & (d.start < r.end) & (d.end > r.start)).values
    ins, out = d[inside], d[~inside]
    if len(ins) and len(out):
        print(f"\n2. PHASING (median phase_pvalue):")
        print(f"   inside PHAS loci : {ins.phase_pvalue.median():.2e}  (n={len(ins):,})")
        print(f"   elsewhere        : {out.phase_pvalue.median():.2e}  (n={len(out):,})")
        print(f"   p<=0.01 inside {(ins.phase_pvalue <= 0.01).mean():.1%} "
              f"against {(out.phase_pvalue <= 0.01).mean():.1%} elsewhere")

    # 3. locus-by-locus agreement
    if ours:
        print(f"\n3. AGREEMENT (reference p vs ours, at the {len(ours)} covered loci):")
        print(f"{'id':<12} {'ref_p':>10} {'our_p':>10} {'class':>10} {'review':>7} {'reads':>9}")
        for i, (pid, rp, op, cl, rv, rd) in enumerate(sorted(ours, key=lambda x: x[1])):
            if i < 12:
                print(f"{pid:<12} {rp:>10.1e} {op:>10.1e} {cl:>10} {rv:>7} {rd:>9.0f}")
        if len(ours) > 12:
            print(f"... and {len(ours) - 12} more")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cls = "phasiRNA"
    for a in sys.argv[1:]:
        if a.startswith("--cls"):
            cls = a.split("=", 1)[1] if "=" in a else "phasiRNA"
    if not args:
        sys.exit(__doc__)
    main(args[0], args[1] if len(args) > 1 else DEFAULT_GFF, cls)
