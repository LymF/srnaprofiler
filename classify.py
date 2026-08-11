"""Orchestrator for the `classify` subcommand."""
import os
from collections import defaultdict
from features import extract_by_chrom
from samio import (ensure_sorted_indexed_bam, mapped_chromosomes, total_mapped,
                   prefers_xm_as_nh, weighted_library_size,
                   reference_alignment_counts)
from mapping import available_memory_mb
from decider import RuleScorer, load_config
from reporter import write_tsv, write_gff3, write_run_json
from utils import log


# Resident bytes per alignment, measured: 339 B (Drosophila), 192 B (Arabidopsis).
_BYTES_PER_ALIGNMENT = 340


DEFAULT_PAD = 500


def resolve_pad(cli_pad, config):
    """CLI wins; otherwise the config's clustering.pad; otherwise the default."""
    if cli_pad is not None:
        return cli_pad
    return (config.get("clustering", {}) or {}).get("pad", DEFAULT_PAD)


def class_size_windows(config):
    """Size windows of the enabled classes, used to group neighbouring islands.

    A disabled class must not reach here: these windows decide which islands may
    merge (`clustering._same_size_group`), so leaving one in would change the
    loci themselves. Adding {24,24} moved an Arabidopsis run from 28,754 to
    32,992 loci with no other change.
    """
    out = []
    for spec in RuleScorer.enabled_classes(config).values():
        size = (spec or {}).get("size")
        if size and "lo" in size and "hi" in size:
            out.append((size["lo"], size["hi"]))
    return out


def mincov_to_depth(mincov_rpm, lib_size):
    """Reads-per-million threshold as an absolute depth for this library."""
    if lib_size <= 0:
        return 1
    return max(1, round(mincov_rpm * lib_size / 1e6))


def safe_thread_count(requested, ref_alignments, available_mb, headroom=0.8):
    """Threads that fit in memory: each worker holds one whole reference, so the
    peak is their sum. Caps rather than warns -- overshooting means the OOM
    killer, not a slowdown.
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


def mark_mirna(records, mirna_intervals, min_fraction=0.5):
    """Label a locus miRNA only when the miRNA covers most of its span."""
    for r in records:
        span = r.end - r.start
        if span <= 0:
            continue
        covered = 0
        for (s, e) in mirna_intervals.get(r.chrom, []):
            covered = max(covered, min(r.end, e) - max(r.start, s))
        if covered / span >= min_fraction:
            r.cls = "miRNA"
            r.confidence = 1.0
            r.review = False
            # The class comes from the annotation, not from evidence measured
            # here, so there is nothing to score. Leaving the 0.0 default would
            # read as "no evidence found", which is a different statement.
            r.score = None


def parse_contig_list(value):
    """Comma-separated contig names, as typed on the command line."""
    if not value:
        return []
    return [name.strip() for name in value.split(",") if name.strip()]


def select_contigs(chroms, exclude):
    """Split the references into the ones to analyse and the ones to drop.

    A name that is not in the BAM is an error: a typo would analyse the
    organelle and report a clean run.
    """
    missing = [name for name in exclude if name not in chroms]
    if missing:
        raise ValueError(
            f"contig(s) not in the BAM: {', '.join(missing)}. "
            f"References with reads: {', '.join(chroms)}")
    dropped = [name for name in chroms if name in set(exclude)]
    kept = [name for name in chroms if name not in set(exclude)]
    return kept, dropped


def run_classify(args):
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)
    config = load_config(args.config)
    ext_cfg = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}
    ext_cfg.update(config.get("extractor", {}))
    ext_cfg["size_windows"] = class_size_windows(config)
    ext_cfg["profile_min_reads"] = (config.get("clustering", {}) or {}).get(
        "profile_min_reads", 2)
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

    if not args.genome:
        log("no -g/--genome: the 5' base null falls back to 0.25 (uniform "
            "composition). Pass the genome so s_1u and s_5a are tested against "
            "each locus's own DNA composition", "WARN")

    min_depth = mincov_to_depth(args.mincov, lib_size)
    log(f"--mincov {args.mincov} RPM -> minimum depth {min_depth} "
        f"for a library of {lib_size:.0f} reads", "INFO")

    mirna_intervals = parse_mirna_bed(args.mirna_loci) if args.mirna_loci else {}
    if mirna_intervals:
        log(f"--mirna-loci: {sum(len(v) for v in mirna_intervals.values())} intervals "
            f"will be carved out as loci of their own", "INFO")

    chroms = mapped_chromosomes(bam_path)
    counts = dict(zip(chroms, reference_alignment_counts(bam_path)))
    exclude = parse_contig_list(getattr(args, "organelle_contigs", None))
    if exclude:
        chroms, dropped = select_contigs(chroms, exclude)
        out = sum(counts[name] for name in dropped)
        total = sum(counts.values())
        pct = 100.0 * out / total if total else 0.0
        log(f"--organelle-contigs: dropping {', '.join(dropped)} "
            f"({out:,} of {total:,} alignments, {pct:.2f}%); the library size "
            f"still counts them", "INFO")

    safe = safe_thread_count(threads, [counts[c] for c in chroms],
                             available_memory_mb())
    if safe < threads:
        log(f"Capping {threads} -> {safe} thread(s): each holds one whole "
            f"chromosome in memory and there is not enough free RAM for more "
            f"(this would be an OOM kill, not a slowdown)", "WARN")
        threads = safe
    log(f"Clustering + features across {len(chroms)} reference(s) with {threads} thread(s)...", "INFO")
    records = extract_by_chrom(bam_path, chroms, lib_size, ext_cfg,
                               min_depth, pad, threads, xm_as_nh, mirna_intervals,
                               genome_path=args.genome)
    log(f"{len(records)} clusters", "INFO")

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
