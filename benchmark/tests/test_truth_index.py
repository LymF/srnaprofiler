"""Integrity of benchmark/truth/truth_index.tsv against the files and the manifest."""
import csv
from pathlib import Path

import pytest

from benchmark.manifest import load_manifest
from benchmark.truth import SEED_DIR, guard_seed_use, load_aliases, load_regions

REPO = Path(__file__).resolve().parents[2]
INDEX = REPO / "benchmark" / "truth" / "truth_index.tsv"
MANIFEST = REPO / "benchmark" / "manifest.tsv"
CLASSES = {"miRNA", "piRNA", "siRNA", "phasiRNA", "phasiRNA24", "hcsiRNA"}


def index_rows() -> list[dict[str, str]]:
    with open(INDEX, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_index_header_and_rows():
    with open(INDEX, newline="") as fh:
        header = next(csv.reader(fh, delimiter="\t"))
    assert header == ["path", "cls", "source", "feature_types", "aliases"]
    assert index_rows()


@pytest.mark.parametrize("row", index_rows(), ids=lambda r: r["path"])
def test_row_loads(row):
    path = REPO / row["path"]
    assert path.is_file(), f"missing truth file {row['path']}"
    assert row["cls"] in CLASSES
    assert row["source"]
    types = {t for t in row["feature_types"].split(",") if t} or None
    aliases = load_aliases(REPO / row["aliases"]) if row["aliases"] else None
    regions = load_regions(path, row["cls"], row["source"], feature_types=types, aliases=aliases)
    assert regions, f"{row['path']} gave no regions"
    assert all(r.start < r.end for r in regions)


def test_paths_unique():
    paths = [r["path"] for r in index_rows()]
    assert len(paths) == len(set(paths))


def test_every_manifest_truth_source_is_indexed():
    indexed = {r["path"] for r in index_rows()}
    for lib in load_manifest(MANIFEST):
        for src in lib.truth_sources:
            assert src in indexed, f"{lib.accession}: {src} has no truth_index row"


def test_every_seed_file_is_indexed():
    indexed = {(REPO / r["path"]).resolve() for r in index_rows()}
    seeds = {p.resolve() for p in SEED_DIR.iterdir() if p.suffix in (".bed", ".gff3")}
    assert seeds <= indexed


def test_no_test_library_lists_seed_truth():
    for lib in load_manifest(MANIFEST):
        for src in lib.truth_sources:
            guard_seed_use(lib, REPO / src)  # raises for test role + seed file
