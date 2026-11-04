"""Curated truth regions (BED / GFF3) with contig checks and a seed-use guard."""
from dataclasses import dataclass
from pathlib import Path

from benchmark.manifest import Library

MIRNA_PRECURSOR_TYPES = {"miRNA_primary_transcript", "pre_miRNA"}
SEED_DIR = Path(__file__).resolve().parent / "truth" / "seed"


@dataclass
class Region:
    """A truth locus; start/end are 0-based half-open (BED convention)."""
    chrom: str
    start: int
    end: int
    cls: str
    source: str
    name: str


def load_aliases(path: Path) -> dict[str, str]:
    """Read a two-column TSV (alias, genome_name)."""
    aliases: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        alias, genome_name = line.split("\t")[:2]
        aliases[alias] = genome_name
    return aliases


def _gff3_name(attrs: str) -> str:
    fields = dict(kv.split("=", 1) for kv in attrs.split(";") if "=" in kv)
    return fields.get("ID") or fields.get("Name") or ""


def load_regions(path: Path, cls: str, source: str, *,
                 feature_types: set[str] | None = None,
                 aliases: dict[str, str] | None = None) -> list[Region]:
    """Load a BED or GFF3 file as 0-based half-open regions."""
    path = Path(path)
    is_gff = path.suffix.lower() in (".gff3", ".gff")
    regions: list[Region] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if is_gff:
            if feature_types is not None and f[2] not in feature_types:
                continue
            start, end = int(f[3]) - 1, int(f[4])
            name = _gff3_name(f[8]) if len(f) > 8 else ""
        else:
            start, end = int(f[1]), int(f[2])
            name = f[3] if len(f) > 3 else ""
        chrom = (aliases or {}).get(f[0], f[0])
        regions.append(Region(chrom, start, end, cls, source, name))
    return regions


def _check_names(regions: list[Region], known: set[str], where: Path) -> None:
    missing = sorted({r.chrom for r in regions} - known)
    if missing:
        raise ValueError(f"contigs absent from {where}: {', '.join(missing)}")


def check_contigs(regions: list[Region], fai: Path) -> None:
    """Raise ValueError listing every region contig missing from the .fai."""
    known = {line.split("\t")[0] for line in Path(fai).read_text().splitlines() if line}
    _check_names(regions, known, fai)


def fasta_contigs(fasta: Path) -> set[str]:
    """First word of every header line of a FASTA (what faidx and the pipeline call a contig)."""
    names = set()
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                words = line[1:].split()
                if words:
                    names.add(words[0])
    return names


def check_contigs_fasta(regions: list[Region], fasta: Path) -> None:
    """Like check_contigs, against the first-word headers of a (raw) FASTA."""
    _check_names(regions, fasta_contigs(fasta), fasta)


def guard_seed_use(lib: Library, path: Path) -> None:
    """Seed truth is tuning-only: refuse it for test libraries."""
    if lib.role == "test" and Path(path).resolve().is_relative_to(SEED_DIR):
        raise ValueError(f"{lib.accession}: seed truth {path} may not be used for test species")
