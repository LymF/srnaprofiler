"""Idempotent download of SRA libraries (via iSeq) and genome FASTA files."""
import gzip
import os
import shutil
import subprocess
import urllib.request
from pathlib import Path
from urllib.parse import urlparse

from benchmark.manifest import Library

FASTA_SUFFIXES = (".fa", ".fasta", ".fna")


def _iseq_bin() -> str:
    return os.environ.get("ISEQ_BIN", "iseq")


def fetch_library(lib: Library, dest: Path) -> Path:
    """Download one single-end run as dest/<accession>.fastq.gz (skipped if present)."""
    dest = Path(dest)
    target = dest / f"{lib.accession}.fastq.gz"
    if target.exists():
        return target
    tmp = dest / f".iseq_{lib.accession}"
    shutil.rmtree(tmp, ignore_errors=True)  # discard leftovers of interrupted runs
    tmp.mkdir(parents=True)
    binary = _iseq_bin()
    env = os.environ.copy()
    if os.sep in binary:  # iSeq needs its sibling tools (axel, pigz, sra-tools) on PATH
        env["PATH"] = str(Path(binary).resolve().parent) + os.pathsep + env.get("PATH", "")
    try:
        subprocess.run([binary, "-i", lib.accession, "-g", "-o", str(tmp)], check=True, env=env)
        files = sorted(tmp.rglob("*.fastq.gz"))
        if len(files) > 1:
            names = ", ".join(f.name for f in files)
            raise ValueError(f"{lib.accession}: paired-end run not supported ({names})")
        if not files:
            raise FileNotFoundError(f"iSeq produced no fastq.gz for {lib.accession}")
        if files[0].name.endswith(("_1.fastq.gz", "_2.fastq.gz", "_R1.fastq.gz", "_R2.fastq.gz")):
            raise ValueError(f"{lib.accession}: paired-end run not supported ({files[0].name})")
        shutil.move(str(files[0]), target)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return target


def _download(url: str, target: Path) -> None:
    with urllib.request.urlopen(url) as resp, open(target, "wb") as out:
        shutil.copyfileobj(resp, out)


def fetch_genome(lib: Library, dest: Path) -> Path:
    """Download lib.genome_url into dest and return an uncompressed .fa path."""
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    name = Path(urlparse(lib.genome_url).path).name
    gzipped = name.endswith(".gz")
    if gzipped:
        name = name[:-3]
    stem = name
    for suffix in FASTA_SUFFIXES:
        if name.endswith(suffix):
            stem = name[: -len(suffix)]
            break
    fasta = dest / f"{stem}.fa"
    if fasta.exists():
        return fasta
    raw = dest / f"{stem}.fa.gz" if gzipped else dest / f"{stem}.fa.part"
    _download(lib.genome_url, raw)
    if gzipped:
        with gzip.open(raw, "rb") as src, open(fasta.with_suffix(".fa.part"), "wb") as out:
            shutil.copyfileobj(src, out)
        raw.unlink()
        raw = fasta.with_suffix(".fa.part")
    raw.rename(fasta)
    return fasta
