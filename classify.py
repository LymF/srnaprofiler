"""Orchestrator for the `classify` subcommand."""
import os
from collections import defaultdict
from features import extract_by_chrom
from samio import ensure_sorted_indexed_bam, mapped_chromosomes, total_mapped
from decider import RuleScorer, load_config
from reporter import write_tsv, write_gff3, write_run_json
from utils import log


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

    base = os.path.join(outdir, args.prefix)
    log(f"Ensuring sorted, indexed BAM (threads={threads})...", "INFO")
    bam_path = ensure_sorted_indexed_bam(args.sam, base, threads)

    lib_size = args.norm
    if lib_size is None:
        lib_size = total_mapped(bam_path)
        log(f"--norm not given: using total mapped reads ({lib_size}) as lib_size for RPM", "WARN")

    chroms = mapped_chromosomes(bam_path)
    log(f"Clustering + features across {len(chroms)} reference(s) with {threads} thread(s)...", "INFO")
    records = extract_by_chrom(bam_path, chroms, lib_size, ext_cfg, args.mincov, args.pad, threads)
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
