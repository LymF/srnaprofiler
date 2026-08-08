#!/usr/bin/env python3
"""Does the score separate real phasing from background of the same size?

Positives: the 8 experimentally characterised TAS loci, and the loci on PHAS21
annotated by sRNAminer. Controls: 21-nt loci on neither, matched to each
positive by distinct_seqs -- unmatched controls would be won by locus size
alone, since the dominant register's mass fraction is 0.954 in loci with <=20
distinct sequences against 1/21 expected.

There is no threshold here that counts as passing. The number is for comparing
formulations: the p-value scored 0.560 (TAS) and 0.481 (PHAS), i.e. chance.

Usage:  python3 extras/validate_phase_evidence.py run.clusters.tsv
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TAS_BED = os.path.join(HERE, "phas_experimental_at_TAIR10.bed")
PHAS_GFF = os.path.join(HERE, "phas21_ref_at_TAIR10.gff3")

# The PHAS annotation is in TAIR names, the runs are on the RefSeq assembly.
TAIR_TO_REFSEQ = {"Chr1": "NC_003070.9", "Chr2": "NC_003071.7",
                  "Chr3": "NC_003074.8", "Chr4": "NC_003075.7",
                  "Chr5": "NC_003076.8"}


def load_intervals():
    tas = pd.read_csv(TAS_BED, sep="\t", comment="#",
                      names=["chrom", "start", "end", "name"])
    phas = []
    with open(PHAS_GFF) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                phas.append((TAIR_TO_REFSEQ.get(f[0], f[0]),
                             int(f[3]) - 1, int(f[4])))
    return tas, pd.DataFrame(phas, columns=["chrom", "start", "end"])


def overlaps(d, iv):
    hit = np.zeros(len(d), dtype=bool)
    for chrom, grp in iv.groupby("chrom"):
        same = (d.chrom == chrom).values
        if not same.any():
            continue
        for _, r in grp.iterrows():
            hit |= same & (d.end.values > r.start) & (d.start.values < r.end)
    return hit


def auc(pos, neg):
    """Probability a positive outranks a control, ties counted as half."""
    a = np.asarray(pos, dtype=float)
    b = np.asarray(neg, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    wins = (a[:, None] > b[None, :]).sum()
    ties = (a[:, None] == b[None, :]).sum()
    return (wins + 0.5 * ties) / (len(a) * len(b))


def main(tsv_path):
    d = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    d = d[d.dicer_call.astype(str) == "21"].reset_index(drop=True)
    tas, phas = load_intervals()
    d["in_tas"] = overlaps(d, tas)
    d["in_phas"] = overlaps(d, phas) | d.in_tas

    pool = d[~d.in_phas]
    negs = []
    for i in d.index[d.in_phas]:
        n = d.at[i, "distinct_seqs"]
        cand = pool[(pool.distinct_seqs >= n * 0.5)
                    & (pool.distinct_seqs <= n * 1.5)]
        cand = cand[~cand.index.isin(negs)]
        if len(cand):
            negs.append(int(cand.sample(1, random_state=len(negs)).index[0]))

    ctrl = d.loc[negs]
    print(f"{len(d)} loci 21-nt; {d.in_tas.sum()} on experimental TAS, "
          f"{d.in_phas.sum()} on annotated PHAS; {len(ctrl)} matched controls")
    print(f"\n{'column':14} {'AUC TAS':>8} {'AUC PHAS':>9}")
    for col in ("score", "s_phase", "phase_number", "phase_ratio"):
        if col not in d.columns:
            continue
        print(f"{col:14} {auc(d[d.in_tas][col], ctrl[col]):8.3f} "
              f"{auc(d[d.in_phas][col], ctrl[col]):9.3f}")
    print("\nmedian score  TAS / PHAS / control: "
          f"{d[d.in_tas].score.median():.3f} / "
          f"{d[d.in_phas].score.median():.3f} / {ctrl.score.median():.3f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1])
