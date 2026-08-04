"""Orchestrator for the `classify` subcommand."""
import os
from collections import defaultdict
from features import extract_by_chrom
from samio import (ensure_sorted_indexed_bam, mapped_chromosomes, total_mapped,
                   prefers_xm_as_nh, weighted_library_size,
                   reference_alignment_counts, mapped_genome_size)
from mapping import available_memory_mb
from decider import RuleScorer, load_config
from reporter import write_tsv, write_gff3, write_run_json
from utils import log


# Resident bytes per alignment, measured: 339 B (Drosophila), 192 B (Arabidopsis).
_BYTES_PER_ALIGNMENT = 340


DEFAULT_PAD = 200


def resolve_pad(cli_pad, config):
    """CLI wins; otherwise the config's clustering.pad; otherwise the default."""
    if cli_pad is not None:
        return cli_pad
    return (config.get("clustering", {}) or {}).get("pad", DEFAULT_PAD)


def mincov_to_depth(mincov_rpm, lib_size):
    """Reads-per-million threshold as an absolute depth for this library."""
    if lib_size <= 0:
        return 1
    return max(1, round(mincov_rpm * lib_size / 1e6))


def safe_thread_count(requested, ref_alignments, available_mb, headroom=0.8):
    """Threads that fit in memory, given one whole reference per worker.

    `extract_by_chrom` uses a multiprocessing Pool, so N workers hold N whole
    references at once and the peak is their SUM. Sizing by the largest
    reference alone would be far too pessimistic on a genome with one big
    chromosome and several small ones, so this adds up the N largest and keeps
    the biggest N that fits.

    It caps rather than warns because overshooting is not a slowdown: the OOM
    killer takes the process out, losing the whole run.
    """
    if isinstance(ref_alignments, int):        # a single count: treat as one reference
        ref_alignments = [ref_alignments]
    counts = sorted((c for c in ref_alignments if c > 0), reverse=True)
    if not counts or available_mb <= 0:
        return max(1, requested)
    budget_mb = available_mb * headroom
    fits, total = 0, 0.0
    for c in counts[:requested]:
        total += (c * _BYTES_PER_ALIGNMENT) / (1024 ** 2)
        if total > budget_mb:
            break
        fits += 1
    return max(1, min(requested, fits))


def parse_mirna_bed(path):
    intervals = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "track")) or not line.strip():
                continue
            f = line.split("\t")
            intervals[f[0]].append((int(f[1]), int(f[2])))
    return intervals


def mark_mirna(records, mirna_intervals):
    for r in records:
        for (s, e) in mirna_intervals.get(r.chrom, []):
            if r.start < e and s < r.end:
                r.cls = "miRNA"
                r.confidence = 1.0
                r.review = False
                break


def run_classify(args):
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    config = load_config(args.config)
    ext_cfg = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}
    ext_cfg.update(config.get("extractor", {}))
    threads = getattr(args, "threads", 1) or 1
    pad = resolve_pad(getattr(args, "pad", None), config)

    base = os.path.join(outdir, args.prefix)
    log(f"Ensuring sorted, indexed BAM (threads={threads})...", "INFO")
    bam_path = ensure_sorted_indexed_bam(args.sam, base, threads)

    xm_as_nh = prefers_xm_as_nh(bam_path)
    if xm_as_nh:
        log("bowtie1 input detected: using XM tag as NH for 1/NH multi-mapper weighting", "INFO")

    lib_size = args.norm
    if lib_size is None:
        # Same 1/NH scale as the per-locus counts, so RPM is not deflated by the
        # alignment lines a multi-mapping read contributes (~4x with bowtie1 -k5).
        log("--norm not given: computing the weighted library size (sum of 1/NH)...", "INFO")
        lib_size = weighted_library_size(bam_path, xm_as_nh)
        log(f"Weighted library size: {lib_size:.1f} reads "
            f"({total_mapped(bam_path)} alignment lines)", "INFO")

    min_depth = mincov_to_depth(args.mincov, lib_size)
    log(f"--mincov {args.mincov} RPM -> minimum depth {min_depth} "
        f"for a library of {lib_size:.0f} reads", "INFO")

    chroms = mapped_chromosomes(bam_path)
    safe = safe_thread_count(threads, reference_alignment_counts(bam_path),
                             available_memory_mb())
    if safe < threads:
        log(f"Capping {threads} -> {safe} thread(s): each holds one whole "
            f"chromosome in memory and there is not enough free RAM for more "
            f"(this would be an OOM kill, not a slowdown)", "WARN")
        threads = safe
    log(f"Clustering + features across {len(chroms)} reference(s) with {threads} thread(s)...", "INFO")
    records = extract_by_chrom(bam_path, chroms, lib_size, ext_cfg,
                               min_depth, pad, threads, xm_as_nh,
                               mapped_genome_size(bam_path))
    log(f"{len(records)} clusters", "INFO")

    mirna_intervals = parse_mirna_bed(args.mirna_loci) if args.mirna_loci else {}
    if mirna_intervals:
        mark_mirna(records, mirna_intervals)

    scorer = RuleScorer(config)
    for r in records:
        if r.cls != "miRNA":
            scorer.decide(r)

    write_tsv(records, base + ".clusters.tsv")
    write_gff3(records, base + ".clusters.gff3")
    params = {k: v for k, v in vars(args).items() if not callable(v)}
    write_run_json(params, base + ".run.json")
    log(f"Written: {base}.clusters.tsv / .gff3 / .run.json", "INFO")
