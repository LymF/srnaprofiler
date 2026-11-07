"""Run one benchmark library end to end: fetch, pipeline, truth, metrics."""
import csv
import json
import shutil
from pathlib import Path

import pysam

from benchmark.fetch import fetch_genome, fetch_library
from benchmark.manifest import Config, Library, load_config, load_manifest
from benchmark.metrics_class import classification_metrics, read_fragments, read_loci
from benchmark.metrics_delim import delimitation_metrics, locus_counts
from benchmark.metrics_expr import expressed_regions
from benchmark.run import result_paths, run_library
from benchmark.truth import (Region, check_contigs, check_contigs_fasta, guard_seed_use,
                             load_aliases, load_regions)

REPO = Path(__file__).resolve().parent.parent
TRUTH_INDEX = REPO / "benchmark" / "truth" / "truth_index.tsv"
CONFIG = REPO / "benchmark" / "config.yaml"
RULES = {"arthropod": "rules.yaml", "plant": "rules.plant.yaml"}


def _index_rows() -> dict[str, dict]:
    with open(TRUTH_INDEX, newline="") as fh:
        return {row["path"]: row for row in csv.DictReader(fh, delimiter="\t")}


def load_truth(lib: Library) -> list[Region]:
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


_load_truth = load_truth


def _library(acc: str, manifest: Path) -> Library:
    libs = {l.accession: l for l in load_manifest(manifest)}
    if acc not in libs:
        raise ValueError(f"unknown accession {acc} in {manifest}")
    return libs[acc]


def _write_metrics(lib: Library, regions: list[Region], cfg: Config, workdir: Path,
                   results_dir: Path) -> Path:
    """Score <workdir>/<ACC>.clusters.tsv against the truth and write metrics.json."""
    acc = lib.accession
    paths = result_paths(workdir, acc)
    std_genome = workdir / f"{acc}.genome.std.fa"
    fai = Path(str(std_genome) + ".fai")
    if not fai.exists():
        pysam.faidx(str(std_genome))
    check_contigs(regions, fai)  # ValueError names the missing contigs

    loci = read_loci(paths.clusters_tsv)
    expressed = expressed_regions(paths.bam, regions, cfg.expressed_min_reads)
    clusters = [r for r in regions if r.cls != "miRNA"]
    sources = sorted({r.source for r in regions})
    cluster_sources = sorted({r.source for r in clusters})
    metrics = {
        "accession": acc,
        "organism": lib.organism,
        "role": lib.role,
        "classification": classification_metrics(loci, regions, cfg, expressed),
        # one delimitation per source: sources overlap each other, so pooling double counts
        "delimitation": {
            "global": locus_counts(loci, cfg, read_fragments(paths.clusters_tsv)),
            "by_source": {
                s: delimitation_metrics(loci, [r for r in clusters if r.source == s], cfg)
                for s in cluster_sources
            },
        },
        "by_source": {
            s: classification_metrics(loci, [r for r in regions if r.source == s], cfg,
                                   expressed)["per_class"]
            for s in sources
        },
        "n_loci": len(loci),
    }
    out = results_dir / acc
    out.mkdir(parents=True, exist_ok=True)
    metrics_path = out / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2))
    shutil.copyfile(paths.meta, out / "benchmark_run.json")
    return metrics_path


def run_benchmark_library(acc: str, manifest: Path, data_dir: Path, results_dir: Path,
                          threads: int, classifier: str = "rules",
                          ml_model: Path | None = None) -> Path:
    """Run one library and return the path of its metrics.json."""
    lib = _library(acc, manifest)
    data_dir, results_dir = Path(data_dir), Path(results_dir)
    cfg = load_config(CONFIG)

    species = lib.organism.replace(" ", "_")
    genome = fetch_genome(lib, data_dir / "genomes" / species)
    # fail on a truth/genome contig mismatch before the long reads download and pipeline
    regions = _load_truth(lib)
    check_contigs_fasta(regions, genome)
    reads = fetch_library(lib, data_dir / "reads")
    rules = REPO / RULES[lib.taxon_group]
    workdir = results_dir / acc / "work"
    run_library(lib, reads, genome, workdir, rules, threads=threads,
                index_prefix=results_dir / "_index" / species / species,
                classifier=classifier, ml_model=ml_model)
    return _write_metrics(lib, regions, cfg, workdir, results_dir)


def rebuild_metrics(acc: str, manifest: Path, results_dir: Path) -> Path:
    """Rebuild metrics.json from an existing results/<ACC>/work, without running anything."""
    lib = _library(acc, manifest)
    results_dir = Path(results_dir)
    workdir = results_dir / acc / "work"
    paths = result_paths(workdir, acc)
    for needed in (paths.clusters_tsv, paths.meta, workdir / f"{acc}.genome.std.fa"):
        if not needed.exists():
            raise FileNotFoundError(f"{acc}: {needed} missing; run the library first")
    return _write_metrics(lib, _load_truth(lib), load_config(CONFIG), workdir, results_dir)
