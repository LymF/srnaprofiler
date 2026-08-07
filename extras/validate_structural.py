#!/usr/bin/env python3
"""Specificity against structural RNA: what must NOT be called a small RNA.

Every other validator asks about recall -- do we find the loci a reference says
are there. None asks the opposite question, which is just as easy to fail:
tRNA, rRNA, snoRNA and snRNA fragments are abundant in any small RNA library
(degradation, and tRNA halves that are real but are not siRNA or piRNA), and a
classifier that labels them siRNA or piRNA is wrong in a way recall cannot see.

There is no threshold here that counts as "passing" -- some overlap is expected
because these genes sit in transcribed neighbourhoods. What the numbers are for
is comparison: between runs, between organisms, and against the share of the
genome the features occupy.

Usage:  python3 extras/validate_structural.py run.clusters.tsv annotation.gff3
"""
import sys
from collections import defaultdict

import pandas as pd

KINDS = ("tRNA", "rRNA", "snoRNA", "snRNA")


def load_structural(gff_path):
    """Structural RNA genes, by kind."""
    rows = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in KINDS:
                continue
            rows.append((f[2], f[0], int(f[3]) - 1, int(f[4])))
    return pd.DataFrame(rows, columns=["kind", "chrom", "start", "end"])


def main(tsv_path, gff_path):
    feat = load_structural(gff_path)
    if feat.empty:
        sys.exit(f"No {'/'.join(KINDS)} feature in {gff_path}")
    d = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    d["span"] = d.end - d.start

    print(f"{len(feat)} structural RNA genes "
          f"({dict(feat.kind.value_counts())}), {len(d):,} loci in the run\n")

    by_chrom = defaultdict(list)
    for r in feat.itertuples():
        by_chrom[r.chrom].append((r.start, r.end, r.kind))

    hit_kind = []
    for r in d.itertuples():
        best = None
        for (s, e, kind) in by_chrom.get(r.chrom, ()):
            if r.start < e and r.end > s:
                cov = min(r.end, e) - max(r.start, s)
                if best is None or cov > best[0]:
                    best = (cov, kind)
        hit_kind.append(best[1] if best else "")
    d = d.assign(structural=hit_kind)

    on = d[d.structural != ""]
    print(f"1. LOCI OVERLAPPING A STRUCTURAL RNA: {len(on):,} "
          f"({len(on) / len(d):.2%} of loci, {on.reads.sum() / d.reads.sum():.2%} "
          f"of the read mass)\n")

    if len(on):
        print("2. WHAT THOSE LOCI ARE CALLED:")
        for cls, n in on.cls.value_counts().items():
            share = n / len(on)
            mass = on.reads[on.cls == cls].sum() / on.reads.sum()
            rev = (on.review[on.cls == cls].astype(str).str.lower() == "true").mean()
            print(f"   {str(cls):<12} {n:>6} loci ({share:5.1%}), {mass:5.1%} of their "
                  f"read mass, {rev:5.1%} flagged review")

        print("\n3. BY FEATURE KIND (loci overlapping each):")
        for kind, grp in on.groupby("structural"):
            top = grp.cls.value_counts()
            first = f"{top.index[0]} {top.iloc[0]}" if len(top) else "-"
            print(f"   {kind:<8} {len(grp):>6} loci, most called {first}")

    print("\n4. FOR COMPARISON — share of each class that sits on a structural RNA:")
    for cls, grp in d.groupby("cls"):
        print(f"   {str(cls):<12} {(grp.structural != '').mean():6.2%} of {len(grp):>7,} loci")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
