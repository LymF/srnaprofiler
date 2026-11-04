"""Truth files added in Task 8c: they load, match PROVENANCE counts, use dm6 contigs and 0-based BED coordinates."""
import csv
import re
from pathlib import Path

import pytest

from benchmark.manifest import load_manifest
from benchmark.truth import SEED_DIR, check_contigs_fasta, guard_seed_use, load_aliases, load_regions

REPO = Path(__file__).resolve().parents[2]
DMEL = REPO / "benchmark" / "truth" / "dmel"
INDEX = REPO / "benchmark" / "truth" / "truth_index.tsv"
# Chromosome arms of Ensembl BDGP6.54 (= dm6); verified against the genome FASTA with check_contigs_fasta (PROVENANCE).
DM6_ARMS = {"2L", "2R", "3L", "3R", "4", "X"}
MIRNA_TYPES = {"miRNA_primary_transcript", "pre_miRNA"}


def _rows() -> dict[str, dict[str, str]]:
    with open(INDEX, newline="") as fh:
        return {r["path"]: r for r in csv.DictReader(fh, delimiter="\t")}


def _load(rel: str):
    row = _rows()[rel]
    types = {t for t in row["feature_types"].split(",") if t} or None
    aliases = load_aliases(REPO / row["aliases"]) if row["aliases"] else None
    return row, load_regions(REPO / rel, row["cls"], row["source"], feature_types=types, aliases=aliases)


def _provenance_count(name: str) -> int:
    """Count stated in PROVENANCE as `name` ... `N loci` on the line(s) after the file heading."""
    text = (REPO / "benchmark" / "PROVENANCE.md").read_text()
    section = text.split(f"`benchmark/truth/dmel/{name}`", 1)[1]
    return int(re.search(r"\*\*(\d+) loci\*\*", section).group(1))


def test_dme_mirbase_loads_with_aliases():
    row, regs = _load("benchmark/truth/dmel/dme_mirbase_v23.1.gff3")
    assert row["cls"] == "miRNA" and row["source"] == "miRBase v23.1"
    assert len(regs) == _provenance_count("dme_mirbase_v23.1.gff3") == 258
    assert {r.chrom for r in regs} <= DM6_ARMS
    assert all(r.cls == "miRNA" for r in regs)


def test_dme_mirbase_raw_names_need_the_alias_file():
    raw = load_regions(DMEL / "dme_mirbase_v23.1.gff3", "miRNA", "s", feature_types=MIRNA_TYPES)
    assert {r.chrom for r in raw}.isdisjoint(DM6_ARMS)
    assert set(load_aliases(DMEL / "dme_contig_aliases.tsv").values()) == DM6_ARMS


def test_dme_endosirna_bed_loads():
    row, regs = _load("benchmark/truth/dmel/dme_endosiRNA_2008_dm6.bed")
    assert row["cls"] == "siRNA"
    assert len(regs) == _provenance_count("dme_endosiRNA_2008_dm6.bed") == 76
    assert {r.chrom for r in regs} <= DM6_ARMS
    assert all(r.cls == "siRNA" and r.name.startswith("Czech2008_S") for r in regs)


def test_dme_endosirna_bed_is_zero_based_half_open():
    for line in (DMEL / "dme_endosiRNA_2008_dm6.bed").read_text().splitlines():
        f = line.split("\t")
        start, end = int(f[1]), int(f[2])
        assert 0 <= start < end, line


def test_dme_endosirna_excludes_mirna_and_piRNA_truth():
    _, endo = _load("benchmark/truth/dmel/dme_endosiRNA_2008_dm6.bed")
    others = []
    for rel in ("benchmark/truth/dmel/dme_mirgenedb_v3.0.gff3", "benchmark/truth/dmel/dme_mirbase_v23.1.gff3",
                "benchmark/truth/dmel/dme_piRNAclusterDB2_proTRAC.gff3"):
        others += _load(rel)[1]
    for e in endo:
        assert not any(o.chrom == e.chrom and o.start < e.end and e.start < o.end for o in others), e


def test_contig_check_helper_accepts_dm6_names(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{c} dna:chromosome\nACGT\n" for c in sorted(DM6_ARMS)))
    _, regs = _load("benchmark/truth/dmel/dme_endosiRNA_2008_dm6.bed")
    check_contigs_fasta(regs, fa)


def test_new_files_are_not_seed_and_no_test_library_uses_seed():
    for rel in ("benchmark/truth/dmel/dme_mirbase_v23.1.gff3", "benchmark/truth/dmel/dme_endosiRNA_2008_dm6.bed"):
        assert not (REPO / rel).resolve().is_relative_to(SEED_DIR)
    for lib in load_manifest(REPO / "benchmark" / "manifest.tsv"):
        for src in lib.truth_sources:
            guard_seed_use(lib, REPO / src)
