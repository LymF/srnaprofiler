"""Run the sRNAProfiler pipeline (qc, mapping, mirna, classify) on one library."""
import hashlib
import json
import os
import shutil
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


def file_md5(path: Path) -> str:
    """md5 of a file, streamed."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pipeline_sha256(repo: Path = REPO) -> str:
    """sha256 over the pipeline *.py files at the repo root (not benchmark/) and rules*.yaml."""
    repo = Path(repo)
    h = hashlib.sha256()
    for f in sorted([*repo.glob("*.py"), *repo.glob("rules*.yaml")], key=lambda f: f.name):
        h.update(f.name.encode() + b"\0" + f.read_bytes() + b"\0")
    return h.hexdigest()


def _git_dirty() -> bool | None:
    """True if tracked files differ from HEAD; None when git is unavailable."""
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO), "status", "--porcelain", "--untracked-files=no"],
            check=True, capture_output=True, text=True,
        )
        return bool(out.stdout.strip())
    except (subprocess.CalledProcessError, OSError):
        return None


def result_paths(workdir: Path, acc: str) -> RunResult:
    """Where a finished run leaves its clusters table, BAM and metadata."""
    workdir = Path(workdir)
    return RunResult(workdir / f"{acc}.clusters.tsv", workdir / f"{acc}.bowtie.bam",
                     workdir / "benchmark_run.json")


def _fresh(outputs: list[Path], inputs: list[Path]) -> bool:
    """All outputs exist and none is older than any input."""
    if not all(o.exists() for o in outputs):
        return False
    newest_input = max((i.stat().st_mtime for i in inputs if i.exists()), default=0.0)
    return min(o.stat().st_mtime for o in outputs) >= newest_input


def _index_ready(prefix: Path) -> bool:
    """A bowtie1 index is exactly six non-empty files; a stray .part file does not count."""
    names = ("1", "2", "3", "4", "rev.1", "rev.2")
    for n in names:
        f = prefix.parent / f"{prefix.name}.{n}.ebwt"
        if not f.is_file() or f.stat().st_size == 0:
            return False
    return True


def _share_index(built: Path, shared: Path) -> None:
    """Copy a freshly built index to the shared prefix, file by file, atomically."""
    shared.parent.mkdir(parents=True, exist_ok=True)
    for f in sorted(built.parent.glob(built.name + ".*.ebwt*")):
        dest = shared.parent / (shared.name + f.name[len(built.name):])
        tmp = dest.with_name(dest.name + ".part")
        shutil.copyfile(f, tmp)
        os.replace(tmp, dest)


def _resolve_ml_model(rules: Path, ml_model: Path | None) -> Path:
    """The model file classify will load: --ml-model, else rules `ml.model` against the repo."""
    if ml_model is not None:
        model = Path(ml_model)
    else:
        import yaml
        cfg = yaml.safe_load(Path(rules).read_text()) or {}
        rel = (cfg.get("ml") or {}).get("model")
        if not rel:
            raise ValueError(f"{rules} has no ml.model; pass --ml-model")
        model = Path(rel)
        if not model.is_absolute():
            model = REPO / model
    if not model.is_file():
        raise FileNotFoundError(f"ML model not found: {model}")
    return model


def run_library(
    lib: Library, reads: Path, genome: Path, workdir: Path, rules: Path,
    threads: int = 8, index_prefix: Path | None = None,
    classifier: str = "rules", ml_model: Path | None = None,
) -> RunResult:
    """Run qc, mapping, mirna and classify; skip a step whose outputs are up to date.

    index_prefix: shared bowtie index prefix for the species. When its index exists mapping
    uses it (--index-prefix) instead of rebuilding; otherwise mapping builds the index in the
    workdir and it is copied there for the next library of the species. mirna uses the same
    prefix; without index_prefix the per-library index is used.
    """
    if classifier not in ("rules", "ml"):
        raise ValueError(f"classifier must be 'rules' or 'ml', got {classifier!r}")
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    model = _resolve_ml_model(Path(rules), ml_model) if classifier == "ml" else None
    reads, genome, rules = Path(reads), Path(genome), Path(rules)
    acc = lib.accession
    py = os.environ.get("SRNAPROFILER_PYTHON", sys.executable)
    base = [py, str(PIPELINE)]
    out = str(workdir)
    t = ["--threads", str(threads)]

    trimmed = workdir / f"{acc}.trimmed.R1.fastq.gz"
    std_genome = workdir / f"{acc}.genome.std.fa"
    local_index = workdir / f"{acc}.bowtie_index"  # built by `mapping -g`
    bam = workdir / f"{acc}.bowtie.bam"
    bai = workdir / f"{acc}.bowtie.bam.bai"
    bed = workdir / f"{acc}.mirna_loci.bed"
    mirna_done = workdir / f"{acc}.mirna.done"  # written after a successful mirna run
    clusters = workdir / f"{acc}.clusters.tsv"
    engine = ENGINES[lib.taxon_group]
    shared = Path(index_prefix) if index_prefix else None

    meta = workdir / "benchmark_run.json"
    previous: dict[str, list[str]] = {}  # commands that produced the outputs we may skip
    previous_sha: str | None = None
    if meta.exists():
        try:
            old = json.loads(meta.read_text())
            previous_sha = old.get("pipeline_sha256")
            previous = {c[2]: c for c in old.get("commands", [])}
        except (ValueError, IndexError, TypeError, AttributeError):
            previous, previous_sha = {}, None
    current_sha = pipeline_sha256()
    # outputs may be reused only if they were made by exactly this pipeline code
    can_resume = previous_sha is not None and previous_sha == current_sha
    commands: dict[str, list[str]] = {}
    skipped: list[str] = []

    def step(name: str, outputs: list[Path], inputs: list[Path], argv: list[str]) -> None:
        if can_resume and _fresh(outputs, inputs):
            skipped.append(name)
            if name in previous:
                commands[name] = previous[name]
        else:
            commands[name] = argv
            subprocess.run(argv, check=True)

    use_shared = shared is not None and _index_ready(shared)
    qc_index = ["--index-prefix", str(shared)] if use_shared else []
    step("qc", [trimmed, std_genome], [reads, genome],
         base + ["qc", "-p", acc, "-o", out, "--reads1", str(reads), "-g", str(genome)]
         + qc_index + t)

    ref = ["--index-prefix", str(shared)] if use_shared else ["-g", str(std_genome)]
    step("mapping", [bam, bai], [trimmed, std_genome],
         base + ["mapping", "-p", acc, "-o", out, "--reads", str(trimmed)] + ref + t)
    if shared is not None and not use_shared and "mapping" not in skipped:
        _share_index(local_index, shared)
    index = shared if shared is not None and _index_ready(shared) else local_index

    # the marker, not the bed, says mirna ran: a library with no miRNA writes no bed
    mirna_out = mirna_done if mirna_done.exists() else bed
    step("mirna", [mirna_out], [trimmed, std_genome],
         base + ["mirna", "-p", acc, "-o", out, "--reads", str(trimmed),
                 "--species", lib.organism, "--engine", engine,
                 "-g", str(std_genome), "--index-prefix", str(index)] + t)
    if "mirna" not in skipped:
        mirna_done.write_text("ok\n")

    use_bed = bed.exists() and bed.stat().st_size > 0
    classify = base + ["classify", "-p", acc, "-o", out, "--sam", str(bam),
                       "-g", str(std_genome), "--config", str(rules)] + t
    if use_bed:
        classify += ["--mirna-loci", str(bed)]
    if model is not None:
        classify += ["--classifier", "ml"]
        if ml_model is not None:
            classify += ["--ml-model", str(ml_model)]
    step("classify", [clusters],
         [bam, rules] + ([bed] if use_bed else []) + ([model] if model else []), classify)

    if not clusters.exists():
        raise FileNotFoundError(f"classify did not produce {clusters}")

    meta.write_text(json.dumps({
        "tool_git_sha": _git_sha(),
        "pipeline_sha256": current_sha,
        "resumed_from_pipeline_sha256": previous_sha,
        "git_dirty": _git_dirty(),
        "reads_md5": file_md5(reads),
        "genome_md5": file_md5(genome),
        "rules_sha256": hashlib.sha256(rules.read_bytes()).hexdigest(),
        "rules_path": str(rules),
        "classifier": classifier,
        "ml_model_path": str(model) if model else None,
        "ml_model_sha256": hashlib.sha256(model.read_bytes()).hexdigest() if model else None,
        "accession": acc,
        "organism": lib.organism,
        "role": lib.role,
        "mirna_engine": engine,
        "mirna_loci_used": use_bed,
        "ncrna_db_passed": False,
        "skipped_steps": skipped,
        "commands": list(commands.values()),
    }, indent=2))
    return RunResult(clusters_tsv=clusters, bam=bam, meta=meta)
