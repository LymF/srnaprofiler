#!/usr/bin/env python3
"""Compare two runs over the same reference loci: a mutant against its wild type.

Recall against an annotation says whether we find what someone else annotated.
A mutant with a known phenotype asks a different question -- whether the number
we report moves when the biology moves. Annotation cannot answer that, because
the annotation is the same file in both runs.

Everything is compared per million weighted reads, on the same reference
intervals, because the confounder here is depth: two libraries never have the
same size, and an unnormalised count would report sequencing effort as biology.

Usage:
    python3 extras/compare_runs.py wt.clusters.tsv mut.clusters.tsv \\
        reference.gff3 FEATURE CLASS [control.gff3 FEATURE CLASS]

The optional second reference is the internal control: a class that should NOT
move. If both move together, the cause is depth or processing, not the pathway.
"""
import sys

import numpy as np
import pandas as pd


def load_intervals(path, kind):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5 and f[2] == kind:
                rows.append((f[0], int(f[3]) - 1, int(f[4])))
    return pd.DataFrame(rows, columns=["chrom", "start", "end"])


def overlaps(d, iv):
    hit = np.zeros(len(d), dtype=bool)
    for chrom, grp in iv.groupby("chrom"):
        same = (d.chrom == chrom).values
        if not same.any():
            continue
        for _, r in grp.iterrows():
            hit |= same & (d.end.values > r.start) & (d.start.values < r.end)
    return hit


def lib_size(d):
    """Weighted library size, recovered from the run's own RPM column."""
    ok = d.reads_rpm > 0
    if not ok.any():
        return float("nan")
    return float(np.median(d.reads[ok] / d.reads_rpm[ok] * 1e6))


def summarise(path, iv, cls):
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d["inside"] = overlaps(d, iv)
    review = d.review.astype(str).str.lower() == "true"
    called = d[(d.cls == cls) & ~review]
    return {
        "lib": lib_size(d),
        "rpm_inside": float(d.loc[d.inside, "reads_rpm"].sum()),
        "loci": int(len(called)),
        "loci_on_ref": int(called["inside"].sum()),
        "rpm_called": float(called.reads_rpm.sum()),
    }


def report(label, wt, mut):
    print(f"\n{label}")
    print(f"  {'':28}{'WT':>14}{'mutante':>14}{'razão':>10}")
    rows = [
        ("massa nos loci de ref. (RPM)", "rpm_inside"),
        (f"loci chamados (não-review)", "loci"),
        ("desses, em loci de ref.", "loci_on_ref"),
        ("massa dos chamados (RPM)", "rpm_called"),
    ]
    for name, key in rows:
        a, b = wt[key], mut[key]
        ratio = (b / a) if a else float("nan")
        fmt = "{:14.0f}" if key.startswith("rpm") else "{:14d}"
        print(f"  {name:28}" + fmt.format(a) + fmt.format(b) + f"{ratio:10.2f}")


def size_profile(path):
    """Share of the library in each Dicer size, from the loci table."""
    d = pd.read_csv(path, sep="\t", low_memory=False)
    out = {}
    for L in ("21", "24"):
        sub = d[d.dicer_call.astype(str) == L]
        out[L] = float(sub.reads_rpm.sum())
    return out


def main(wt_path, mut_path, ref, ref_kind, cls, *control):
    iv = load_intervals(ref, ref_kind)
    if iv.empty:
        sys.exit(f"no {ref_kind} feature in {ref}")
    wt, mut = summarise(wt_path, iv, cls), summarise(mut_path, iv, cls)
    print(f"WT : {wt_path}\n     library {wt['lib']:,.0f} weighted reads")
    print(f"MUT: {mut_path}\n     library {mut['lib']:,.0f} weighted reads")
    report(f"{cls} against {ref_kind} ({len(iv)} reference loci)", wt, mut)

    if len(control) >= 3:
        c_iv = load_intervals(control[0], control[1])
        c_wt = summarise(wt_path, c_iv, control[2])
        c_mut = summarise(mut_path, c_iv, control[2])
        report(f"CONTROLE -- {control[2]} against {control[1]} "
               f"({len(c_iv)} loci), which should NOT move", c_wt, c_mut)

    print("\nmassa por tamanho dominante (RPM somado sobre os loci):")
    a, b = size_profile(wt_path), size_profile(mut_path)
    for L in ("21", "24"):
        ratio = b[L] / a[L] if a[L] else float("nan")
        print(f"  {L} nt{'':22}{a[L]:14.0f}{b[L]:14.0f}{ratio:10.2f}")


if __name__ == "__main__":
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
