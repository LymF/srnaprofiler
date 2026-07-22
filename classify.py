"""Orquestrador do subcomando `classify`."""
import os
from collections import defaultdict
from clustering import cluster_regions
from features import extract
from samio import read_alignments
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

    log("Definindo clusters por cobertura...", "INFO")
    regions = cluster_regions(args.sam, args.mincov, args.pad)
    log(f"{len(regions)} clusters", "INFO")

    lib_size = args.norm
    if lib_size is None:
        lib_size = sum(1 for _ in read_alignments(args.sam))
        log(f"--norm ausente: usando total de reads mapeados ({lib_size}) como lib_size para RPM", "WARN")

    records = extract(args.sam, regions, lib_size, ext_cfg)

    mirna_intervals = parse_mirna_bed(args.mirna_loci) if args.mirna_loci else {}
    if mirna_intervals:
        mark_mirna(records, mirna_intervals)

    scorer = RuleScorer(config)
    for r in records:
        if r.cls != "miRNA":
            scorer.decide(r)

    base = os.path.join(outdir, args.prefix)
    write_tsv(records, base + ".clusters.tsv")
    write_gff3(records, base + ".clusters.gff3")
    write_run_json(vars(args), base + ".run.json")
    log(f"Escrito: {base}.clusters.tsv / .gff3 / .run.json", "INFO")
