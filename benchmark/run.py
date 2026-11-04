"""Run the sRNAProfiler pipeline (qc, mapping, mirna, classify) on one library."""
import hashlib
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from benchmark.manifest import Library

REPO = Path(__file__).resolve().parent.parent
PIPELINE = REPO / "srnaprofiler.py"
ENGINES = {"arthropod": "mirdeep2", "plant": "mirdeep-p2"}


@dataclass
class RunResult:
    clusters_tsv: Path
    bam: Path
    meta: Path


def _git_sha() -> str:
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        )
        return out.stdout.strip() or "unknown"
    except (subprocess.CalledProcessError, OSError):
        return "unknown"


def run_library(
    lib: Library, reads: Path, genome: Path, workdir: Path, rules: Path,
    threads: int = 8,
) -> RunResult:
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    acc = lib.accession
    py = os.environ.get("SRNAPROFILER_PYTHON", sys.executable)
    base = [py, str(PIPELINE)]
    out = str(workdir)
    t = ["--threads", str(threads)]

    trimmed = workdir / f"{acc}.trimmed.R1.fastq.gz"
    std_genome = workdir / f"{acc}.genome.std.fa"
    index = workdir / f"{acc}.bowtie_index"  # built by `mapping -g`
    bam = workdir / f"{acc}.bowtie.bam"
    bed = workdir / f"{acc}.mirna_loci.bed"
    clusters = workdir / f"{acc}.clusters.tsv"
    engine = ENGINES[lib.taxon_group]

    commands: list[list[str]] = []

    def run(argv: list[str]) -> None:
        commands.append(argv)
        subprocess.run(argv, check=True)

    run(base + ["qc", "-p", acc, "-o", out, "--reads1", str(reads),
                "-g", str(genome)] + t)
    run(base + ["mapping", "-p", acc, "-o", out, "--reads", str(trimmed),
                "-g", str(std_genome)] + t)
    run(base + ["mirna", "-p", acc, "-o", out, "--reads", str(trimmed),
                "--species", lib.organism, "--engine", engine,
                "-g", str(std_genome), "--index-prefix", str(index)] + t)

    use_bed = bed.exists() and bed.stat().st_size > 0
    classify = base + ["classify", "-p", acc, "-o", out, "--sam", str(bam),
                       "-g", str(std_genome), "--config", str(rules)] + t
    if use_bed:
        classify += ["--mirna-loci", str(bed)]
    run(classify)

    if not clusters.exists():
        raise FileNotFoundError(f"classify did not produce {clusters}")

    meta = workdir / "benchmark_run.json"
    meta.write_text(json.dumps({
        "tool_git_sha": _git_sha(),
        "rules_sha256": hashlib.sha256(Path(rules).read_bytes()).hexdigest(),
        "rules_path": str(rules),
        "accession": acc,
        "organism": lib.organism,
        "role": lib.role,
        "mirna_engine": engine,
        "mirna_loci_used": use_bed,
        "ncrna_db_passed": False,
        "commands": commands,
    }, indent=2))
    return RunResult(clusters_tsv=clusters, bam=bam, meta=meta)
