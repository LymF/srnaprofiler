#!/usr/bin/env python3
"""Validate a run against the curated miRNAs of miRBase.

The other validators only ask about one class each -- piRNA clusters in
Drosophila, PHAS loci in plants -- and both references are themselves
predictions. miRBase is curated, and it covers both kingdoms, so it is the
strongest reference available for a class the tool must not get wrong.

It answers two questions, and the second is the one nobody was asking:

  1. recall      -- annotated miRNA precursors covered by a locus of ours
  2. SPECIFICITY -- what those loci are called when `--mirna-loci` is NOT given.
                    A miRNA hairpin turning up as piRNA or siRNA is a false
                    positive of that class, and recall-only validation is blind
                    to it.

miRBase GFF3 coordinates match the assemblies already in use (dme on Release_6,
ath on TAIR10), so no liftover is involved:

    curl -sSLO https://www.mirbase.org/download/dme.gff3

Usage:  python3 extras/validate_mirna.py run.clusters.tsv dme.gff3
"""
import re
import sys

import pandas as pd

# Chromosome naming differs per source; map to a comparable key.
REFSEQ = {
    # Drosophila Release 6
    "NT_033779": "2l", "NT_033778": "2r", "NC_004354": "x",
    "NT_037436": "3l", "NT_033777": "3r", "NC_004353": "4",
    # Arabidopsis TAIR10
    "NC_003070": "1", "NC_003071": "2", "NC_003074": "3",
    "NC_003075": "4", "NC_003076": "5",
}


def norm_chrom(c):
    s = str(c)
    base = s.split(".")[0]
    if base in REFSEQ:
        return REFSEQ[base]
    return re.sub(r"^chr", "", s.lower())


def load_mirna(gff_path):
    """miRNA precursors (the hairpin, not the mature arm)."""
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "miRNA_primary_transcript":
                continue
            name = dict(re.findall(r"(\w+)=([^;]+)", f[8])).get("Name", "")
            rows.append((name, f[0], int(f[3]) - 1, int(f[4])))
    return pd.DataFrame(rows, columns=["name", "chrom", "start", "end"])


def main(tsv_path, gff_path):
    mir = load_mirna(gff_path)
    if mir.empty:
        sys.exit(f"No miRNA_primary_transcript in {gff_path}")
    d = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    d = d.assign(k=d.chrom.map(norm_chrom))
    mir = mir.assign(k=mir.chrom.map(norm_chrom))
    if not set(mir.k) & set(d.k):
        sys.exit(f"No chromosome in common. miRBase: {sorted(set(mir.k))[:6]} "
                 f"run: {sorted(set(d.k))[:6]}")

    print(f"{len(mir)} annotated miRNA precursors, {len(d):,} loci in the run\n")

    covered, calls, review_n = 0, [], 0
    for r in mir.itertuples():
        same = d[d.k == r.k]
        ov = same[(same.start < r.end) & (same.end > r.start)]
        if not len(ov):
            continue
        covered += 1
        best = ov.loc[ov.reads.idxmax()]
        calls.append(best.cls)
        review_n += str(best.review).lower() == "true"

    print(f"1. RECALL: {covered}/{len(mir)} precursors are covered by a locus "
          f"({covered / len(mir):.1%})")
    if not calls:
        return
    counts = pd.Series(calls).value_counts()
    print(f"\n2. WHAT THEY ARE CALLED (no --mirna-loci; {review_n} of {covered} "
          f"are flagged for review):")
    for cls, n in counts.items():
        label = "correct" if cls == "miRNA" else "FALSE POSITIVE of this class"
        print(f"   {str(cls):<12} {n:>5}  ({n / covered:5.1%})  {label}")

    # How much of each class is explained by miRNA hairpins? The other side of
    # the same question: a class polluted by miRNA is a class we cannot trust.
    inside = pd.Series(False, index=d.index)
    for r in mir.itertuples():
        inside |= ((d.k == r.k) & (d.start < r.end) & (d.end > r.start))
    print("\n3. SHARE OF EACH CLASS SITTING ON AN ANNOTATED miRNA:")
    for cls, grp in d.groupby("cls"):
        share = inside[grp.index].mean()
        print(f"   {str(cls):<12} {share:6.2%} of {len(grp):>7,} loci")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
