import dataclasses
from pathlib import Path

import pytest

from benchmark.manifest import Library, load_config, load_manifest, validate_coverage

REAL_CONFIG = Path(__file__).resolve().parents[1] / "config.yaml"

DM = "Drosophila melanogaster"
AT = "Arabidopsis thaliana"
OS = "Oryza sativa"


def lib(accession, organism, role, taxon_group, **kw):
    defaults = dict(tissue="whole", replicate=1, genome_url="https://example.org/g.fa.gz", truth_sources=[])
    defaults.update(kw)
    return Library(accession, organism, role, taxon_group, **defaults)


def good_set():
    return [
        lib("SRR1", DM, "tuning", "arthropod"),
        lib("SRR2", AT, "tuning", "plant"),
        lib("SRR3", OS, "tuning", "plant"),
        lib("SRR4", "Bombyx mori", "test", "arthropod"),
        lib("SRR5", "Zea mays", "test", "plant"),
    ]


def write_tsv(path, libs):
    header = [f.name for f in dataclasses.fields(Library)]
    lines = ["\t".join(header)]
    for l in libs:
        row = []
        for name in header:
            v = getattr(l, name)
            row.append(",".join(v) if isinstance(v, list) else str(v))
        lines.append("\t".join(row))
    path.write_text("\n".join(lines) + "\n")


def test_config_defaults_match_spec():
    cfg = load_config(REAL_CONFIG)
    assert cfg.tiny_bp == [50, 100]
    assert cfg.low_reads == [20, 50, 100]
    assert cfg.min_overlap_frac == 0.5


def test_load_manifest_roundtrip(tmp_path):
    libs = good_set()
    libs[0] = lib("SRR1", DM, "tuning", "arthropod", replicate=2, truth_sources=["a.bed", "b.gff3"])
    p = tmp_path / "manifest.tsv"
    write_tsv(p, libs)
    assert load_manifest(p) == libs


@pytest.mark.parametrize("col,value", [("role", "holdout"), ("taxon_group", "fungus")])
def test_load_manifest_rejects_unknown_values(tmp_path, col, value):
    p = tmp_path / "manifest.tsv"
    write_tsv(p, good_set())
    lines = p.read_text().splitlines()
    idx = lines[0].split("\t").index(col)
    cells = lines[1].split("\t")
    cells[idx] = value
    lines[1] = "\t".join(cells)
    p.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match=col):
        load_manifest(p)


def test_validate_coverage_accepts_good_set():
    validate_coverage(good_set())


def test_validate_coverage_rejects_missing_test_species():
    no_arthropod = [l for l in good_set() if l.accession != "SRR4"]
    with pytest.raises(ValueError, match="arthropod"):
        validate_coverage(no_arthropod)
    no_plant = [l for l in good_set() if l.accession != "SRR5"]
    with pytest.raises(ValueError, match="plant"):
        validate_coverage(no_plant)
    # a tuning species in the test role does not count as a distinct test species
    only_dm = no_arthropod + [lib("SRR6", DM, "test", "arthropod")]
    with pytest.raises(ValueError):
        validate_coverage(only_dm)


@pytest.mark.parametrize("organism,group", [(OS, "plant"), (AT, "plant"), (DM, "arthropod")])
def test_rice_cannot_be_test_role(organism, group):
    libs = good_set() + [lib("SRR9", organism, "test", group)]
    with pytest.raises(ValueError, match=organism):
        validate_coverage(libs)


@pytest.mark.parametrize("missing", [DM, AT, OS])
def test_tuning_species_required(missing):
    libs = [l for l in good_set() if l.organism != missing]
    with pytest.raises(ValueError, match=missing):
        validate_coverage(libs)


def test_check_manifest_cli(tmp_path, capsys):
    from benchmark.__main__ import main

    p = tmp_path / "manifest.tsv"
    write_tsv(p, good_set())
    assert main(["check-manifest", "--manifest", str(p)]) == 0
    assert "OK: 5 libraries, coverage satisfied" in capsys.readouterr().out

    write_tsv(p, good_set()[:3])
    assert main(["check-manifest", "--manifest", str(p)]) == 1
    assert capsys.readouterr().err.strip()
