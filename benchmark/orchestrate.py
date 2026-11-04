"""Run one benchmark library end to end: fetch, pipeline, truth, metrics."""
import csv
import json
import shutil
from pathlib import Path

import pysam

from benchmark.fetch import fetch_genome, fetch_library
from benchmark.manifest import Library, load_config, load_manifest
from benchmark.metrics_class import classification_metrics, read_loci
from benchmark.metrics_delim import delimitation_metrics, locus_counts
from benchmark.run import run_library
from benchmark.truth import (Region, check_contigs, check_contigs_fasta, guard_seed_use,
                             load_aliases, load_regions)

REPO = Path(__file__).resolve().parent.parent
TRUTH_INDEX = REPO / "benchmark" / "truth" / "truth_index.tsv"
CONFIG = REPO / "benchmark" / "config.yaml"
RULES = {"arthropod": "rules.yaml", "plant": "rules.plant.yaml"}


def _index_rows() -> dict[str, dict]:
    with open(TRUTH_INDEX, newline="") as fh:
        return {row["path"]: row for row in csv.DictReader(fh, delimiter="\t")}


def _load_truth(lib: Library) -> list[Region]:
    index = _index_rows()
    regions: list[Region] = []
    for path in lib.truth_sources:
        row = index.get(path)
        if row is None:
            raise ValueError(f"{lib.accession}: truth source {path} has no row in {TRUTH_INDEX}")
        full = REPO / path
        guard_seed_use(lib, full)
        types = {t for t in row["feature_types"].split(",") if t} or None
        aliases = load_aliases(REPO / row["aliases"]) if row["aliases"] else None
        regions += load_regions(full, row["cls"], row["source"],
                                feature_types=types, aliases=aliases)
    return regions


def run_benchmark_library(acc: str, manifest: Path, data_dir: Path, results_dir: Path,
                          threads: int) -> Path:
    """Run one library and return the path of its metrics.json."""
    libs = {l.accession: l for l in load_manifest(manifest)}
    if acc not in libs:
        raise ValueError(f"unknown accession {acc} in {manifest}")
    lib = libs[acc]
    data_dir, results_dir = Path(data_dir), Path(results_dir)
    cfg = load_config(CONFIG)

    genome = fetch_genome(lib, data_dir / "genomes" / lib.organism.replace(" ", "_"))
    # fail on a truth/genome contig mismatch before the long reads download and pipeline
    regions = _load_truth(lib)
    check_contigs_fasta(regions, genome)
    reads = fetch_library(lib, data_dir / "reads")
    rules = REPO / RULES[lib.taxon_group]
    workdir = results_dir / acc / "work"
    res = run_library(lib, reads, genome, workdir, rules, threads=threads)

    std_genome = workdir / f"{acc}.genome.std.fa"
    fai = Path(str(std_genome) + ".fai")
    if not fai.exists():
        pysam.faidx(str(std_genome))
    check_contigs(regions, fai)  # ValueError names the missing contigs

    loci = read_loci(res.clusters_tsv)
    clusters = [r for r in regions if r.cls != "miRNA"]
    sources = sorted({r.source for r in regions})
    cluster_sources = sorted({r.source for r in clusters})
    metrics = {
        "accession": acc,
        "organism": lib.organism,
        "role": lib.role,
        "classification": classification_metrics(loci, regions, cfg),
        # one delimitation per source: sources overlap each other, so pooling double counts
        "delimitation": {
            "global": locus_counts(loci, cfg),
            "by_source": {
                s: delimitation_metrics(loci, [r for r in clusters if r.source == s], cfg)
                for s in cluster_sources
            },
        },
        "by_source": {
            s: classification_metrics(loci, [r for r in regions if r.source == s], cfg)["per_class"]
            for s in sources
        },
        "n_loci": len(loci),
    }
    out = results_dir / acc
    out.mkdir(parents=True, exist_ok=True)
    metrics_path = out / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2))
    shutil.copyfile(res.meta, out / "benchmark_run.json")
    return metrics_path
