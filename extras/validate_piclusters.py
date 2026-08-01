#!/usr/bin/env python3
"""Validate classify calls against known D. melanogaster piRNA clusters.

Cross-checks a `.clusters.tsv` against `extras/piC_ref_dm6.bed` (conserved ovary
piRNA clusters in Release 6 coordinates) and reports:

  1. recall     -- reference clusters holding at least one non-review piRNA call
  2. enrichment -- piRNA vs siRNA read mass inside the reference clusters,
                   against what the clusters occupy of the genome
  3. landmarks  -- flamenco against 42AB, the discriminative test: a uni-strand
                   somatic cluster and a dual-strand germline one have opposite
                   published signatures, and nothing in the pipeline knows where
                   either of them is

Reference coordinates: Srivastav, Feschotte & Clark (2024), Genome Res
34(5):711-724, https://doi.org/10.1101/gr.278062.123

Usage:  python3 extras/validate_piclusters.py run.clusters.tsv [piC_ref_dm6.bed]
"""
import os
import sys
import numpy as np
import pandas as pd

# RefSeq names of the Release 6 chromosome arms (the BED uses UCSC-style names).
CHROM = {"2L": "NT_033779.5", "2R": "NT_033778.4", "3L": "NT_037436.4",
         "3R": "NT_033777.3", "X": "NC_004354.4", "4": "NC_004353.4"}
ARM_LEN = {"NT_033779.5": 23513712, "NT_033778.4": 25286936, "NT_037436.4": 28110227,
           "NT_033777.3": 32079331, "NC_004354.4": 23542271, "NC_004353.4": 1348131}

# Two landmarks with opposite, well-established signatures.
LANDMARKS = {
    "flamenco (X, uni-strand/somatic)": ("NC_004354.4", 21624500, 21995076),
    "42AB (2R, dual-strand/germline)": ("NT_033778.4", 6236220, 6515294),
}


def load_reference(path):
    ref = pd.read_csv(path, sep="\t", comment="#",
                      names=["arm", "start", "end", "n_strains"])
    ref["chrom"] = ref.arm.map(CHROM)
    return ref.dropna(subset=["chrom"])


def mark_inside(clusters, ref):
    """Boolean per locus: does it overlap any reference cluster?"""
    inside = np.zeros(len(clusters), bool)
    for chrom, g in ref.groupby("chrom"):
        idx = np.where((clusters.chrom == chrom).values)[0]
        if not len(idx):
            continue
        s, e = clusters.start.values[idx], clusters.end.values[idx]
        hit = np.zeros(len(idx), bool)
        for rs, re_ in g[["start", "end"]].values:
            hit |= (s < re_) & (e > rs)
        inside[idx] = hit
    return inside


def main(tsv_path, bed_path):
    ref = load_reference(bed_path)
    d = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    d = d[d.chrom.isin(ARM_LEN)]
    if not len(d):
        sys.exit("No locus on a main chromosome arm — is this a Release 6 run?")
    review = d.review.astype(str).str.lower() == "true"
    d = d.assign(inside=mark_inside(d, ref))

    ref_bp = sum(min(e, ARM_LEN[c]) - s for c, s, e in ref[["chrom", "start", "end"]].values)
    genome_bp = sum(ARM_LEN.values())
    frac_genome = ref_bp / genome_bp
    print(f"reference: {len(ref)} piCs, {ref_bp / 1e6:.2f} Mb "
          f"= {frac_genome * 100:.2f}% of the main arms\n")

    pirna = d[(d.cls == "piRNA") & ~review.values]
    recovered = 0
    for chrom, g in ref.groupby("chrom"):
        p = pirna[pirna.chrom == chrom]
        for rs, re_ in g[["start", "end"]].values:
            if ((p.start < re_) & (p.end > rs)).any():
                recovered += 1
    print(f"1. RECALL: {recovered}/{len(ref)} reference piCs "
          f"({recovered / len(ref) * 100:.1f}%) hold >=1 non-review piRNA locus")

    print(f"\n2. ENRICHMENT (read mass inside the piCs, against {frac_genome * 100:.2f}% "
          f"of the genome):")
    for cls in ("piRNA", "siRNA"):
        m = d.cls == cls
        if not m.any():
            continue
        share = d.reads[m & d.inside].sum() / d.reads[m].sum()
        print(f"   {cls}: {share * 100:5.1f}%  -> {share / frac_genome:.1f}x")

    print("\n3. LANDMARKS (read-weighted, expected: flamenco uni-strand and "
          "ping-pong-poor, 42AB dual-strand and ping-pong-rich):")
    for name, (c, s, e) in LANDMARKS.items():
        g = d[(d.chrom == c) & (d.start < e) & (d.end > s)]
        if not len(g):
            print(f"   {name}: no locus")
            continue
        w = g.reads / g.reads.sum()
        print(f"   {name}")
        print(f"      loci={len(g):4d}  reads={g.reads.sum():10.0f}  "
              f"piRNA={(g.cls == 'piRNA').mean():.3f}")
        print(f"      frac_top={np.average(g.frac_top, weights=w):.3f}  "
              f"|0.5-frac_top|={np.average((g.frac_top - 0.5).abs(), weights=w):.3f}  "
              f"pingpong_frac={np.average(g.pingpong_frac, weights=w):.3f}  "
              f"1U={np.average(g.pct_1u, weights=w):.3f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    default_bed = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "piC_ref_dm6.bed")
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else default_bed)
