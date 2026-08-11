"""Every deliverable, named the same way whichever engine ran."""
from pathlib import Path
from types import SimpleNamespace

from mirna import finish_mirna, deliverable_paths

PREDS = [
    {"chrom": "Chr1", "strand": "+", "name": "p1", "start": 1234, "end": 1325,
     "mature": "TGGACTGAAGGGAGCTCCCTT"},
    {"chrom": "Chr2", "strand": "-", "name": "p2", "start": 88100, "end": 88190,
     "mature": "TTGACAGAAGATAGAGAGCAC"},
]


def args_for(tmp_path):
    return SimpleNamespace(prefix="at", outdir=str(tmp_path), engine="mirdeep-p2",
                           genome=None, no_validate=True, threads=1)


def test_the_deliverables_are_named_off_the_prefix(tmp_path):
    p = deliverable_paths(tmp_path, "at")
    assert p["table"].name == "at.mirna.tsv"
    assert p["bed"].name == "at.mirna_loci.bed"
    assert p["mature"].name == "at.mirna.mature.fa"
    assert p["star"].name == "at.mirna.star.fa"
    assert p["precursor"].name == "at.mirna.precursor.fa"
    assert p["run_json"].name == "at.mirna.run.json"
    assert p["structures"].name == "at.mirna_structures"


def test_without_validation_the_table_is_written_and_no_fasta_is(tmp_path):
    """--no-validate: every row NotValidated, so there is no validated set."""
    finish_mirna(PREDS, args_for(tmp_path), tmp_path / "work", reads=None)
    p = deliverable_paths(tmp_path, "at")
    lines = p["table"].read_text().splitlines()
    assert len(lines) == 3
    assert all("NotValidated" in ln for ln in lines[1:])
    assert not p["mature"].exists()
    assert not p["star"].exists()
    assert not p["precursor"].exists()


def test_the_bed_is_still_written_unvalidated_so_the_run_can_go_on(tmp_path):
    bed = finish_mirna(PREDS, args_for(tmp_path), tmp_path / "work", reads=None)
    assert bed == deliverable_paths(tmp_path, "at")["bed"]
    assert len(bed.read_text().splitlines()) == 2


def test_the_run_json_records_what_ran(tmp_path):
    import json
    finish_mirna(PREDS, args_for(tmp_path), tmp_path / "work", reads=None)
    meta = json.loads(deliverable_paths(tmp_path, "at")["run_json"].read_text())
    assert meta["engine"] == "mirdeep-p2"
    assert meta["validated"] is False
    assert meta["counts"]["predicted"] == 2
    assert meta["counts"]["passed"] == 0


def validated_args(tmp_path):
    return SimpleNamespace(prefix="at", outdir=str(tmp_path), engine="mirdeep-p2",
                           genome="genome.fa", no_validate=False, threads=1)


def fake_validation(monkeypatch, verdict):
    """miRScore's answer for p1, without miRScore: it lives in another env."""
    from mirna import Validation
    rows = {"p1": {"result": verdict, "msSeq": "GUUCC", "precSeq": "UUCUUGG",
                   "flags": "NA"}}
    monkeypatch.setattr(
        "mirna.validate_with_mirscore",
        lambda *a, **k: Validation(rows=rows, submitted={"p1", "p2"}))


def test_a_rerun_that_validates_nothing_leaves_no_fasta_from_the_old_one(
        tmp_path, monkeypatch):
    """The four files describe one set, so an empty FASTA beats a stale one."""
    fake_validation(monkeypatch, "Pass")
    finish_mirna(PREDS, validated_args(tmp_path), tmp_path / "work", reads=None)
    p = deliverable_paths(tmp_path, "at")
    assert ">p1" in p["mature"].read_text()

    fake_validation(monkeypatch, "Fail")
    finish_mirna(PREDS, validated_args(tmp_path), tmp_path / "work", reads=None)
    assert p["bed"].read_text() == ""
    for key in ("mature", "star", "precursor"):
        assert p[key].read_text() == "", f"{key} still claims the old pass set"


def test_nothing_submitted_is_NA_and_not_a_measured_zero(tmp_path):
    import json
    finish_mirna(PREDS, args_for(tmp_path), tmp_path / "work", reads=None)
    meta = json.loads(deliverable_paths(tmp_path, "at")["run_json"].read_text())
    assert meta["counts"]["submitted"] == "NA"
