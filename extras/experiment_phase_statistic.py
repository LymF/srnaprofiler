#!/usr/bin/env python3
"""Does any abundance-aware phasing statistic beat the current one?

Truth sets (Arabidopsis): the 8 experimentally characterised TAS loci and the
loci overlapping annotated PHAS21. Negatives: 21-nt loci far from any PHAS
annotation, matched to the positives by distinct_seqs -- unmatched controls
would be won by locus size alone, which is the trap the raw phase_ratio falls
into (median 0.954 among loci with <=20 distinct sequences, against 1/21).

Candidates, all on the same window the current test picks:
  hyper        current: hypergeometric over positions, Bonferroni for the scan
  ratio        raw mass fraction in the dominant register (what sRNAminer reports)
  perm_ratio   the same fraction, but as a permutation p-value
  perm_score   permutation p of phase_score (ratio x cycles x ln(mass))
  perm_number  permutation p of the number of cycles in the dominant register

The permutation keeps the abundances and the number of reads and throws the
positions uniformly over the window, so it prices in the selection of the
maximum over 21 registers -- which is what no test here or in the field does.
"""
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

sys.path.insert(0, "/home/lucas/srnaprofiler")
from samio import read_region_handle, prefers_xm_as_nh
from signals import phase_test, phase_metrics

BAM = "/mnt/d/srnaprofiler-teste/plant/classify/at_ncrna.sorted.bam"
TSV = "/mnt/d/srnaprofiler-teste/plant/classify/at_ncrna_org.clusters.tsv"
TAS_BED = "/home/lucas/srnaprofiler/extras/phas_experimental_at_TAIR10.bed"
PHAS_GFF = "/home/lucas/srnaprofiler/extras/phas21_ref_at_TAIR10.gff3"
PHASE = 21
CYCLES = 11
B = 2000
RNG = np.random.default_rng(0)

TAIR_TO_REFSEQ = {"Chr1": "NC_003070.9", "Chr2": "NC_003071.7", "Chr3": "NC_003074.8",
                  "Chr4": "NC_003075.7", "Chr5": "NC_003076.8"}


def phase_positions(fh, chrom, start, end, xm, phase_len):
    """5' positions of reads of the phase length, pooled across strands.

    Same construction as features.build_record: the - strand's 5' end is moved
    by -2 so Dicer's 2-nt 3' overhang puts both strands in one register.
    """
    out = []
    for a in read_region_handle(fh, chrom, start, end, xm):
        if len(a.seq) != phase_len:
            continue
        w = 1.0 / a.nh
        if a.strand == "+":
            out.append((a.pos, w))
        else:
            out.append((a.pos + len(a.seq) - 1 - 2, w))
    return out


def best_window(positions, region_start, phase_len):
    """The window the current test would report, as (positions, span_start)."""
    if not positions:
        return [], region_start
    span = max(p for p, _ in positions) - min(p for p, _ in positions) + 1
    wlen = CYCLES * phase_len
    if span <= wlen:
        return positions, region_start
    ordered = sorted(positions)
    coords = [p for p, _ in ordered]
    best, best_key, best_start = None, None, region_start
    import bisect
    for s in range(coords[0], coords[-1] + 1, phase_len):
        i = bisect.bisect_left(coords, s)
        j = bisect.bisect_left(coords, s + wlen)
        if i == j:
            continue
        sub = ordered[i:j]
        score, number, _ = phase_metrics(sub, region_start, phase_len)
        key = (number, score)
        if best_key is None or key > best_key:
            best, best_key, best_start = sub, key, s
    return (best if best is not None else positions), best_start


def observed_stats(win, region_start, phase_len):
    score, number, ratio = phase_metrics(win, region_start, phase_len)
    return dict(score=score, number=number, ratio=ratio)


def permutation_pvalues(win, region_start, phase_len, win_start):
    """Null: same abundances, same read count, positions thrown uniformly."""
    if not win:
        return dict(perm_ratio=1.0, perm_score=1.0, perm_number=1.0)
    by_pos = defaultdict(float)
    for p_, a in win:
        by_pos[p_] += a
    coords = np.array(sorted(by_pos), dtype=np.int64)
    w = np.array([by_pos[c] for c in coords], dtype=float)
    n = len(w)
    obs = observed_stats(win, region_start, phase_len)

    # Keep the occupied positions, shuffle which one carries which mass. The
    # hypergeometric already asks whether the in-register POSITIONS are too many;
    # this asks the complementary question, whether the MASS sits on them.
    reg = np.tile((coords - region_start) % phase_len, (B, 1))
    cyc = np.tile((coords - region_start) // phase_len, (B, 1))
    w = np.array([RNG.permutation(w) for _ in range(B)])

    ratios = np.empty(B)
    numbers = np.empty(B)
    scores = np.empty(B)
    total = w[0].sum()
    for b in range(B):
        mass = np.bincount(reg[b], weights=w[b], minlength=phase_len)
        dom = int(mass.argmax())
        m = mass[dom]
        sel = (reg[b] == dom) & (w[b] > 0)
        cycles = np.unique(cyc[b][sel]).size
        ratios[b] = m / total if total else 0.0
        numbers[b] = cycles
        scores[b] = (m / total) * cycles * np.log(m) if m > 0 else 0.0

    def p(null, value):
        return (1 + np.sum(null >= value)) / (B + 1)

    return dict(perm_ratio=p(ratios, obs["ratio"]),
                perm_score=p(scores, obs["score"]),
                perm_number=p(numbers, obs["number"]))


def load_intervals():
    tas = pd.read_csv(TAS_BED, sep="\t", comment="#",
                      names=["chrom", "start", "end", "name"])
    phas = []
    if os.path.exists(PHAS_GFF):
        for line in open(PHAS_GFF):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                phas.append((TAIR_TO_REFSEQ.get(f[0], f[0]), int(f[3]) - 1, int(f[4])))
    return tas, pd.DataFrame(phas, columns=["chrom", "start", "end"])


def overlaps(d, iv):
    hit = np.zeros(len(d), dtype=bool)
    for chrom, grp in iv.groupby("chrom"):
        m = (d.chrom == chrom).values
        if not m.any():
            continue
        s = d.start.values
        e = d.end.values
        for _, r in grp.iterrows():
            hit |= m & (e > r.start) & (s < r.end)
    return hit


def main():
    d = pd.read_csv(TSV, sep="\t", low_memory=False)
    d = d[d.dicer_call.astype(str) == "21"].reset_index(drop=True)
    tas, phas = load_intervals()
    d["in_tas"] = overlaps(d, tas)
    d["in_phas"] = overlaps(d, phas) | d.in_tas
    print(f"{len(d)} loci 21-nt; {d.in_tas.sum()} sobre TAS experimental, "
          f"{d.in_phas.sum()} sobre PHAS anotado")

    pos_idx = d.index[d.in_phas].tolist()
    # matched negatives: same distinct_seqs decile, never in an annotation
    neg_pool = d[~d.in_phas]
    negs = []
    for i in pos_idx:
        n = d.at[i, "distinct_seqs"]
        cand = neg_pool[(neg_pool.distinct_seqs >= n * 0.5) &
                        (neg_pool.distinct_seqs <= n * 1.5)]
        cand = cand[~cand.index.isin(negs)]
        if len(cand):
            negs.append(int(cand.sample(1, random_state=len(negs)).index[0]))
    print(f"{len(pos_idx)} positivos, {len(negs)} controles casados")

    xm = prefers_xm_as_nh(BAM)
    rows = []
    with pysam.AlignmentFile(BAM, "rb") as fh:
        for label, idxs in (("pos", pos_idx), ("neg", negs)):
            for i in idxs:
                r = d.loc[i]
                p = phase_positions(fh, r.chrom, int(r.start), int(r.end), xm, PHASE)
                if not p:
                    continue
                win, wstart = best_window(p, int(r.start), PHASE)
                stats = observed_stats(win, int(r.start), PHASE)
                perm = permutation_pvalues(win, int(r.start), PHASE, wstart)
                rows.append(dict(label=label, name=str(r["name"]),
                                 is_tas=bool(r.in_tas), reads=r.reads,
                                 distinct=r.distinct_seqs,
                                 hyper=r.phase_pvalue, ratio=stats["ratio"],
                                 number=stats["number"], **perm))
    out = pd.DataFrame(rows)
    out.to_csv("phase_experiment.csv", index=False)
    print(f"\n{len(out)} loci avaliados\n")

    def auc(col, higher_is_better=False):
        a = out[out.label == "pos"][col].values
        b = out[out.label == "neg"][col].values
        if len(a) == 0 or len(b) == 0:
            return float("nan")
        s = 1.0 if higher_is_better else -1.0
        wins = (s * a[:, None] > s * b[None, :]).sum()
        ties = (a[:, None] == b[None, :]).sum()
        return (wins + 0.5 * ties) / (len(a) * len(b))

    print(f"{'estatistica':14} {'AUC':>6}  {'recall pos @ FPR 5%':>20}")
    for col, hib in [("hyper", False), ("ratio", True), ("perm_ratio", False),
                     ("perm_score", False), ("perm_number", False)]:
        a = out[out.label == "pos"][col].values
        b = out[out.label == "neg"][col].values
        if hib:
            thr = np.quantile(b, 0.95)
            rec = (a > thr).mean()
        else:
            thr = np.quantile(b, 0.05)
            rec = (a < thr).mean()
        print(f"{col:14} {auc(col, hib):6.3f}  {rec*100:19.1f}%")

    print("\nnos 8 TAS experimentais (p <= 0.01, ou ratio no top-5% dos controles):")
    t = out[out.is_tas]
    for col, hib in [("hyper", False), ("perm_ratio", False),
                     ("perm_score", False), ("perm_number", False)]:
        print(f"  {col:12} {(t[col] <= 0.01).sum()}/{len(t)}")


if __name__ == "__main__":
    main()
