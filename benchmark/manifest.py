"""Benchmark manifest and config loading, plus species-coverage validation."""
import csv
import dataclasses
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import yaml

ROLES = ("tuning", "test")
TAXON_GROUPS = ("arthropod", "plant")

# Species used for tuning; they must never appear with role "test".
TUNING_SPECIES = {
    "Drosophila melanogaster": "arthropod",
    "Arabidopsis thaliana": "plant",
    "Oryza sativa": "plant",
}


@dataclass
class Library:
    accession: str
    organism: str
    role: Literal["tuning", "test"]
    taxon_group: Literal["arthropod", "plant"]
    tissue: str
    replicate: int
    genome_url: str
    truth_sources: list[str]


@dataclass
class Config:
    tiny_bp: list[int]
    low_reads: list[int]
    min_overlap_frac: float


def load_config(path: Path) -> Config:
    data = yaml.safe_load(Path(path).read_text())
    return Config(
        tiny_bp=[int(x) for x in data["tiny_bp"]],
        low_reads=[int(x) for x in data["low_reads"]],
        min_overlap_frac=float(data["min_overlap_frac"]),
    )


def load_manifest(path: Path) -> list[Library]:
    """Read a tab-separated manifest whose header lists the Library field names."""
    names = [f.name for f in dataclasses.fields(Library)]
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames != names:
            raise ValueError(f"manifest header must be {names}, got {reader.fieldnames}")
        libs = []
        for n, row in enumerate(reader, start=2):
            if row["role"] not in ROLES:
                raise ValueError(f"line {n}: unknown role {row['role']!r}")
            if row["taxon_group"] not in TAXON_GROUPS:
                raise ValueError(f"line {n}: unknown taxon_group {row['taxon_group']!r}")
            try:
                replicate = int(row["replicate"])
            except ValueError:
                raise ValueError(f"line {n}: replicate must be an integer, got {row['replicate']!r}")
            truth = [s.strip() for s in row["truth_sources"].split(",") if s.strip()]
            libs.append(Library(
                accession=row["accession"], organism=row["organism"], role=row["role"],
                taxon_group=row["taxon_group"], tissue=row["tissue"], replicate=replicate,
                genome_url=row["genome_url"], truth_sources=truth,
            ))
    return libs


def validate_coverage(libs: list[Library]) -> None:
    """Raise ValueError unless tuning/test species coverage rules are met."""
    for lib in libs:
        if lib.role == "test" and lib.organism in TUNING_SPECIES:
            raise ValueError(f"{lib.organism} is a tuning species and cannot have role 'test' ({lib.accession})")

    tuning = {l.organism for l in libs if l.role == "tuning"}
    for species in TUNING_SPECIES:
        if species not in tuning:
            raise ValueError(f"missing tuning library for {species}")

    test = [l for l in libs if l.role == "test"]
    for group in TAXON_GROUPS:
        if not any(l.taxon_group == group for l in test):
            raise ValueError(f"test set needs at least one {group} outside the tuning species")
