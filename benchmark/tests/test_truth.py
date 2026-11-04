from pathlib import Path

import pytest

from benchmark.manifest import Library
from benchmark.truth import (MIRNA_PRECURSOR_TYPES, SEED_DIR, Region, check_contigs,
                             guard_seed_use, load_aliases, load_regions)

TRUTH_DIR = Path(__file__).resolve().parents[1] / "truth"


def make_lib(role):
    return Library("SRR1", "Test species", role, "plant", "leaf", 1, "http://x", [])


def test_load_bed_regions(tmp_path):
    p = tmp_path / "a.bed"
    p.write_text("# comment\n\n2L\t100\t200\tclusterA\t0\t+\n3R\t5\t9\n")
    regs = load_regions(p, "piRNA", "test")
    assert regs[0] == Region("2L", 100, 200, "piRNA", "test", "clusterA")
    assert regs[1].name == "" and len(regs) == 2


def test_load_gff3_regions(tmp_path):
    p = tmp_path / "a.gff3"
    p.write_text("##gff-version 3\n"
                 "1\t.\tmiRNA_primary_transcript\t11\t20\t.\t+\t.\tID=MI1;Name=x\n"
                 "1\t.\tmiRNA\t11\t15\t.\t+\t.\tID=MIMAT1\n"
                 "1\t.\tpre_miRNA\t31\t40\t.\t-\t.\tName=only_name\n")
    regs = load_regions(p, "miRNA", "mirbase", feature_types=MIRNA_PRECURSOR_TYPES)
    assert [(r.start, r.end, r.name) for r in regs] == [(10, 20, "MI1"), (30, 40, "only_name")]
    assert len(load_regions(p, "miRNA", "mirbase")) == 3


def test_aliases_rename_contigs(tmp_path):
    p = tmp_path / "a.bed"
    p.write_text("Chr1\t1\t5\tn\n")
    a = tmp_path / "al.tsv"
    a.write_text("Chr1\t1\n")
    assert load_regions(p, "c", "s", aliases=load_aliases(a))[0].chrom == "1"


def test_contig_mismatch_raises(tmp_path):
    fai = tmp_path / "g.fa.fai"
    fai.write_text("2L\t1000\t0\t60\t61\n3R\t1000\t0\t60\t61\n")
    regs = [Region("chr2L", 0, 5, "c", "s", "a"), Region("2L", 0, 5, "c", "s", "b"),
            Region("chrX", 0, 5, "c", "s", "c")]
    with pytest.raises(ValueError, match="chr2L") as e:
        check_contigs(regs, fai)
    assert "chrX" in str(e.value)
    check_contigs(regs[1:2], fai)


def test_seed_files_rejected_for_test_species():
    seed = SEED_DIR / "piC_ref_dm6.bed"
    with pytest.raises(ValueError):
        guard_seed_use(make_lib("test"), seed)
    guard_seed_use(make_lib("tuning"), seed)
    guard_seed_use(make_lib("test"), TRUTH_DIR / "bmor" / "x.gff3")


def test_seed_aliases_map_refseq_and_chr_names():
    a = load_aliases(SEED_DIR / "contig_aliases.tsv")
    assert a["NC_003071.7"] == "2" and a["Chr5"] == "5" and a["NC_003076.8"] == "5"
    regs = load_regions(SEED_DIR / "phas_experimental_at_TAIR10.bed", "PHAS", "seed", aliases=a)
    assert {r.chrom for r in regs} <= {"1", "2", "3", "4", "5"}


@pytest.mark.parametrize("path", sorted(
    p for p in TRUTH_DIR.rglob("*") if p.suffix in (".bed", ".gff3") and SEED_DIR not in p.parents),
    ids=lambda p: p.name)
def test_committed_truth_is_sane(path):
    types = MIRNA_PRECURSOR_TYPES if "mir" in path.name else None  # precursor filter only for miRNA files
    regs = load_regions(path, "x", "y", feature_types=types)
    assert regs
    assert all(r.start < r.end for r in regs)
