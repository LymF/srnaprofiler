"""Extração de features por cluster → ClusterRecord."""
from collections import Counter
from typing import Optional
from records import ClusterRecord
from samio import read_alignments

_COMP = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _region_alns(region, alns):
    chrom, start, end = region
    out = []
    for a in alns:
        if a.chrom == chrom and start <= a.pos < end:
            out.append(a)
    return out


def build_record(region, alns, lib_size, config):
    chrom, start, end = region
    dmin, dmax = config["dicermin"], config["dicermax"]
    scut = config["strand_cutoff"]
    r = ClusterRecord(chrom=chrom, start=start, end=end, name=f"{chrom}:{start}-{end}")

    local = _region_alns(region, alns)
    r.reads = len(local)
    if r.reads == 0:
        return r

    strand_count = Counter()
    seqs = Counter()
    seq_lens = Counter()
    gc_bases = 0
    total_bases = 0
    n_1u = n_10a = 0
    multimap = 0
    phase_ready = []  # (pos5, oriented_len) usado na Task 6

    for a in local:
        oriented = a.seq if a.strand == "+" else revcomp(a.seq)
        L = len(oriented)
        strand_count[a.strand] += 1
        seqs[oriented] += 1
        if L < dmin:
            seq_lens["short"] += 1
        elif L > dmax:
            seq_lens["long"] += 1
        else:
            seq_lens[str(L)] += 1
        gc_bases += sum(1 for b in oriented if b in "GCgc")
        total_bases += L
        if oriented[:1] in ("T", "U"):
            n_1u += 1
        if L >= 10 and oriented[9] == "A":
            n_10a += 1
        if a.nh > 1:
            multimap += 1

    r.distinct_seqs = len(seqs)
    r.major_rna, r.major_rna_reads = seqs.most_common(1)[0]
    r.frac_top = round(strand_count["+"] / r.reads, 4)
    if r.frac_top >= scut:
        r.strand_call = "+"
    elif r.frac_top <= 1 - scut:
        r.strand_call = "-"
    else:
        r.strand_call = "."
    r.short = seq_lens["short"]
    r.long = seq_lens["long"]
    predom_len, predom_reads = None, -1
    for L in range(dmin, dmax + 1):
        c = seq_lens.get(str(L), 0)
        r.size_counts[L] = c
        if c > predom_reads:
            predom_len, predom_reads = L, c
    in_window = r.reads - r.short - r.long
    r.frac_in_window = round(in_window / r.reads, 4)
    r.dicer_call = str(predom_len) if r.frac_in_window >= 0.8 else "N"
    r.pct_1u = round(n_1u / r.reads, 4)
    r.pct_10a = round(n_10a / r.reads, 4)
    r.gc = round(gc_bases / total_bases, 4) if total_bases else 0.0
    r.multimap_frac = round(multimap / r.reads, 4)
    r.reads_rpm = round(r.reads / lib_size * 1e6, 4) if lib_size else 0.0
    return r


def extract(sam_path, regions, lib_size, config):
    alns = list(read_alignments(sam_path))
    return [build_record(region, alns, lib_size, config) for region in regions]
