"""Per-cluster feature extraction → ClusterRecord."""
import bisect
from collections import Counter, defaultdict
from records import ClusterRecord
from samio import read_alignments, read_region, reference_lengths
from clustering import call_islands, merge_compatible, carve_out
from signals import offset_spectrum, pingpong_zscore, phase_test
from evidence import (base_bias_score, strand_bias_score, strand_balance_score,
                      strand_bidirectional_score,
                      pingpong_score_empirical, score_from_p)

_COMP = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def build_record(region, region_alns, lib_size, config):
    """Builds a ClusterRecord from the reads ALREADY assigned to the region.

    `region_alns` must contain all reads that overlap the region by span
    (assignment done in `extract`); build_record does not re-filter by position.
    """
    chrom, start, end = region
    dmin, dmax = config["dicermin"], config["dicermax"]
    scut = config["strand_cutoff"]
    r = ClusterRecord(chrom=chrom, start=start, end=end, name=f"{chrom}:{start}-{end}")

    local = region_alns
    raw_n = len(local)
    if raw_n == 0:
        return r

    # Multi-mapping reads are counted fractionally by 1/NH: a read mapping to N
    # places contributes 1/N at each locus, so its total across the genome is 1.
    # This keeps counts at repetitive loci (e.g. transposon piRNA clusters) honest
    # instead of inflating them once per reported alignment.
    strand_count = Counter()
    seqs = Counter()
    seq_lens = Counter()
    gc_bases = 0.0
    total_bases = 0.0
    n_1u = n_10a = n_5a = 0.0
    multimap = 0
    unique_w = 0.0

    for a in local:
        w = 1.0 / a.nh
        oriented = a.seq if a.strand == "+" else revcomp(a.seq)
        L = len(oriented)
        strand_count[a.strand] += w
        seqs[oriented] += w
        if L < dmin:
            seq_lens["short"] += w
        elif L > dmax:
            seq_lens["long"] += w
        else:
            seq_lens[str(L)] += w
        gc_bases += w * sum(1 for b in oriented if b in "GCgc")
        total_bases += w * L
        if oriented[:1] in ("T", "U"):
            n_1u += w
        elif oriented[:1] == "A":
            # 5'A -> AGO4 in plants (24-nt hc-siRNA). Not pct_10a, which is base 10.
            n_5a += w
        if L >= 10 and oriented[9] == "A":
            n_10a += w
        if a.nh > 1:
            multimap += 1
        else:
            unique_w += w

    r.reads = round(strand_count["+"] + strand_count["-"], 4)
    r.distinct_seqs = len(seqs)
    major_rna, major_w = seqs.most_common(1)[0]
    r.major_rna, r.major_rna_reads = major_rna, round(major_w, 4)
    r.frac_top = round(strand_count["+"] / r.reads, 4)
    if r.frac_top >= scut:
        r.strand_call = "+"
    elif r.frac_top <= 1 - scut:
        r.strand_call = "-"
    else:
        r.strand_call = "."
    r.short = round(seq_lens["short"], 4)
    r.long = round(seq_lens["long"], 4)
    predom_len, predom_reads = None, -1.0
    for L in range(dmin, dmax + 1):
        c = seq_lens.get(str(L), 0)
        r.size_counts[L] = round(c, 4)
        if c > predom_reads:
            predom_len, predom_reads = L, c
    in_window = r.reads - r.short - r.long
    r.frac_in_window = round(in_window / r.reads, 4)
    r.dicer_call = str(predom_len) if r.frac_in_window >= 0.8 else "N"
    in_window_by_len = {int(k): v for k, v in seq_lens.items() if k.isdigit()}
    win_total = sum(in_window_by_len.values())
    if win_total:
        ranked = sorted(in_window_by_len.items(), key=lambda kv: kv[1], reverse=True)
        r.size_purity = round(ranked[0][1] / win_total, 4)
        # A second size is only worth reporting when no single population dominates.
        if len(ranked) > 1 and r.size_purity < 0.8:
            r.dicer_call_2 = str(ranked[1][0])
    r.pct_1u = round(n_1u / r.reads, 4)
    r.pct_10a = round(n_10a / r.reads, 4)
    r.pct_5a = round(n_5a / r.reads, 4)
    # n = distinct_seqs, not the 1/NH read count: PCR duplicates and stacked
    # identical reads are not independent trials, so the read count
    # pseudoreplicates the binomial n. k is scaled to the same denominator.
    r.s_1u = base_bias_score(r.pct_1u * r.distinct_seqs, r.distinct_seqs)
    r.s_5a = base_bias_score(r.pct_5a * r.distinct_seqs, r.distinct_seqs)
    r.s_fractop = strand_bias_score(r.frac_top, r.distinct_seqs)
    r.s_fractop_bal = strand_balance_score(r.frac_top, r.distinct_seqs)
    r.s_fractop_bid = strand_bidirectional_score(r.frac_top, r.distinct_seqs)
    r.gc = round(gc_bases / total_bases, 4) if total_bases else 0.0
    # multimap_frac is a RAW diagnostic: fraction of alignment lines with NH>1.
    r.multimap_frac = round(multimap / raw_n, 4)
    # unique_frac is on the 1/NH scale of `reads`: how much of the locus's mass is
    # evidence of its own, rather than borrowed from repeated copies elsewhere.
    r.unique_frac = round(unique_w / r.reads, 4) if r.reads else 0.0
    r.reads_rpm = round(r.reads / lib_size * 1e6, 4) if lib_size else 0.0

    # --- biogenesis signals (5' positions weighted by 1/NH, as above) ---
    plus5 = Counter()
    minus5 = Counter()
    phase_positions = []
    for a in local:
        w = 1.0 / a.nh
        L = len(a.seq)
        if a.strand == "+":
            five = a.pos
            phase5 = five
        else:
            five = a.pos + L - 1
            # 2-nt 3' overhang correction: Dicer's phased products have a
            # 2-nt 3' overhang, so the 5' end of the - strand sits +2 relative
            # to the + strand's register at the same cut. Subtract 2 to pool both
            # strands into the SAME phase register. (-2 sign fixed by biology/test.)
            phase5 = five - 2
        (plus5 if a.strand == "+" else minus5)[five] += w
        phase_positions.append((phase5, w, L))
    # Wide enough to reach the duplex offset (L-2) of any size in the window;
    # ping-pong keeps its own 1..20 window so its scale does not shift.
    wide_counts = offset_spectrum(plus5, minus5, lo=1, hi=dmax)
    overlap_counts = Counter({d: c for d, c in wide_counts.items() if d <= 20})
    r.pingpong_z, r.pingpong_pairs10 = pingpong_zscore(overlap_counts)
    total_pairs = sum(overlap_counts.values())
    # Dicer duplex: the two strands of one duplex leave their 5' ends L-2 apart.
    wide_total = sum(wide_counts.values())
    if r.dicer_call != "N" and wide_total:
        r.overlap_dup_frac = round(
            wide_counts.get(int(r.dicer_call) - 2, 0) / wide_total, 4)
    # 10-nt 5'-5' overlap fraction: bounded [0,1] and comparable across loci (unlike
    # the raw z). Guarded by a minimum number of pairs so noise can't yield a spurious 1.0.
    r.pingpong_frac = round(overlap_counts.get(10, 0) / total_pairs, 4) if total_pairs >= 10 else 0.0
    r.s_pingpong = pingpong_score_empirical(overlap_counts)

    if r.dicer_call != "N":
        pl = int(r.dicer_call)
        # Only reads of the phase length carry the register (as PhasiHunter).
        in_phase_len = [(p, w) for p, w, L in phase_positions if L == pl]
        res = phase_test(in_phase_len, region_start=start, phase_length=pl)
        r.phase_score, r.phase_number = res.score, res.number
        r.phase_ratio, r.phase_pvalue = res.ratio, res.pvalue
        r.s_phase = score_from_p(r.phase_pvalue)
    return r


def _assign_reads(regions, alns):
    """Assigns each read to the regions it overlaps BY SPAN, in a single pass.

    Overlap by span: `a.pos < region.end and a.pos + len(a.seq) > region.start`
    — consistent with the coverage built in `clustering.cluster_regions`.
    Regions from `cluster_regions` are disjoint and sorted by chromosome, which
    lets us locate the block of overlapping regions via binary search instead
    of rescanning every read per region (avoids O(regions × reads)).
    """
    buckets = [[] for _ in regions]
    per_chrom = defaultdict(list)   # chrom -> [(start, end, orig_index)] sorted by start
    for i, (c, s, e) in enumerate(regions):
        per_chrom[c].append((s, e, i))
    for c in per_chrom:
        per_chrom[c].sort()
    starts = {c: [t[0] for t in rs] for c, rs in per_chrom.items()}

    for a in alns:
        rs = per_chrom.get(a.chrom)
        if not rs:
            continue
        aend = a.pos + len(a.seq)
        st = starts[a.chrom]
        # regions with start < aend are indices [0, j); walk backward
        # while end > a.pos (ends are also increasing → contiguous block).
        j = bisect.bisect_left(st, aend)
        k = j - 1
        while k >= 0:
            s, e, oi = rs[k]
            if e <= a.pos:
                break
            buckets[oi].append(a)
            k -= 1
    return buckets


def extract(sam_path, regions, lib_size, config):
    alns = list(read_alignments(sam_path))
    buckets = _assign_reads(regions, alns)
    return [build_record(region, buckets[i], lib_size, config)
            for i, region in enumerate(regions)]


def chromosome_depth(alns):
    """Per-position depth on the 1/NH scale, matching the per-locus counts.

    The mincov threshold is derived from the weighted library size, so counting
    alignment lines here would make it looser by the average multiplicity.
    """
    depth = defaultdict(float)
    for a in alns:
        w = 1.0 / a.nh
        for p in range(a.pos, a.pos + len(a.seq)):
            depth[p] += w
    return depth


def island_profiles(islands, alns, min_reads=0, size_windows=None,
                    min_window_frac=0.5):
    """Size profile and strand fraction per island, for the merge rule.

    The size call is the class window holding most of the read mass, which
    survives the handful of reads a raw island has -- the modal length does not,
    and flips with one read either way. No window over `min_window_frac` of the
    mass means unknown, and unknown never blocks a merge.
    """
    sizes = [Counter() for _ in islands]
    plus = [0.0] * len(islands)
    total = [0.0] * len(islands)
    starts = [s for s, _ in islands]
    for a in alns:
        j = bisect.bisect_right(starts, a.pos) - 1
        if j < 0 or a.pos >= islands[j][1]:
            continue
        w = 1.0 / a.nh
        sizes[j][len(a.seq)] += w
        total[j] += w
        if a.strand == "+":
            plus[j] += w
    out = []
    for j in range(len(islands)):
        if not total[j] or total[j] < min_reads:
            out.append({"size": None, "window": None, "frac_top": None})
            continue
        window = None
        if size_windows:
            mass = [sum(v for L, v in sizes[j].items() if lo <= L <= hi)
                    for lo, hi in size_windows]
            best = max(range(len(mass)), key=lambda i: mass[i])
            if mass[best] / total[j] > min_window_frac:
                window = best
        out.append({"size": sizes[j].most_common(1)[0][0], "window": window,
                    "frac_top": round(plus[j] / total[j], 4)})
    return out


def process_chromosome(bam_path, chrom, mincov, pad, lib_size, config,
                       xm_as_nh=False, chrom_len=0, carve=()):
    """Cluster and extract features for a single chromosome of an indexed BAM.

    Reads only this chromosome's alignments into memory (via fetch), builds its
    coverage islands, assigns reads by span, and returns one ClusterRecord per
    island. This is the unit of parallelism (chromosomes are independent).
    """
    alns = list(read_region(bam_path, chrom, xm_as_nh=xm_as_nh))
    if not alns:
        return []
    depth = chromosome_depth(alns)
    islands = call_islands(depth, mincov, 0)
    profiles = island_profiles(islands, alns,
                               min_reads=config.get("profile_min_reads", 0),
                               size_windows=config.get("size_windows"))
    merged = merge_compatible(islands, profiles, pad,
                              chrom_len=chrom_len or None,
                              size_windows=config.get("size_windows"))
    # Known miRNA loci become loci of their own instead of being swallowed.
    regions = [(chrom, s, e) for (s, e) in carve_out(merged, carve)]
    if not regions:
        return []
    buckets = _assign_reads(regions, alns)
    return [build_record(regions[i], buckets[i], lib_size, config)
            for i in range(len(regions))]


def _chrom_job(job):
    """Top-level wrapper so multiprocessing can pickle the worker call."""
    return process_chromosome(*job)


def extract_by_chrom(bam_path, chroms, lib_size, config, mincov, pad, threads=1,
                     xm_as_nh=False, carve_by_chrom=None):
    """Cluster + extract features across chromosomes, optionally in parallel.

    Each chromosome is processed independently by `process_chromosome`. With
    `threads > 1` and more than one chromosome, the work is spread over a
    multiprocessing pool; otherwise it runs serially. Records are returned in
    chromosome order.
    """
    lengths = reference_lengths(bam_path)
    carve_by_chrom = carve_by_chrom or {}
    jobs = [(bam_path, c, mincov, pad, lib_size, config, xm_as_nh,
             lengths.get(c, 0), tuple(carve_by_chrom.get(c, ()))) for c in chroms]
    if threads > 1 and len(jobs) > 1:
        import multiprocessing as mp
        with mp.Pool(min(threads, len(jobs))) as pool:
            results = pool.map(_chrom_job, jobs)
    else:
        results = [_chrom_job(j) for j in jobs]
    records = []
    for rs in results:
        records.extend(rs)
    return records
