"""The plots miRScore already makes, moved where someone will find them."""
from pathlib import Path

from mirna_report import build_records, collect_structures

# Real-shaped names: dots and underscores inside, and a longer sibling that a
# careless glob would sweep up with p1.
P1, P2, P3 = "NC_003070.9_142", "NC_003070.9_1420", "NC_003070.9_1484"

PREDS = [
    {"chrom": "NC_003070.9", "strand": "+", "name": P1, "start": 10, "end": 100,
     "mature": "TGGACTGAAG"},
    {"chrom": "NC_003070.9", "strand": "+", "name": P2, "start": 200, "end": 300,
     "mature": "TTGACAGAAG"},
    {"chrom": "NC_003070.9", "strand": "+", "name": P3, "start": 400, "end": 500,
     "mature": "TAAGGGATCT"},
]
ROWS = {P1: {"result": "Pass"}, P2: {"result": "Fail"}}


def recs():
    return build_records(PREDS, "mirdeep-p2", ROWS, submitted={P1, P2})


def make_mirscore_out(tmp_path, subdir, suffix):
    """The two real layouts: strucVis writes <name>.ps.pdf, RNAplots <name>_plot.pdf."""
    out = tmp_path / "mirscore"
    for verdict in ("passed", "failed"):
        (out / subdir / verdict).mkdir(parents=True)
    (out / subdir / "passed" / f"{P1}{suffix}").write_text("plot p1")
    (out / subdir / "failed" / f"{P2}{suffix}").write_text("plot p2")
    return out


def test_strucvis_plots_are_split_by_verdict(tmp_path):
    out = make_mirscore_out(tmp_path, "strucVis", ".ps.pdf")
    dest = tmp_path / "at.mirna_structures"
    assert collect_structures(out, recs(), dest) == 2
    assert (dest / "passed" / f"{P1}.ps.pdf").read_text() == "plot p1"
    assert (dest / "failed" / f"{P2}.ps.pdf").read_text() == "plot p2"


def test_rnaplots_are_the_fallback_when_strucvis_is_absent(tmp_path):
    out = make_mirscore_out(tmp_path, "RNAplots", "_plot.pdf")
    dest = tmp_path / "at.mirna_structures"
    assert collect_structures(out, recs(), dest) == 2
    assert (dest / "passed" / f"{P1}_plot.pdf").exists()


def test_a_longer_name_is_not_swept_up_by_a_shorter_one(tmp_path):
    """NC_003070.9_142 must not collect NC_003070.9_1420's plot."""
    out = make_mirscore_out(tmp_path, "strucVis", ".ps.pdf")
    (out / "strucVis" / "passed" / f"{P2}.ps.pdf").write_text("plot p2 passed")
    dest = tmp_path / "at.mirna_structures"
    collect_structures(out, recs(), dest)
    assert not (dest / "passed" / f"{P2}.ps.pdf").exists()


def test_strucvis_wins_when_both_exist(tmp_path):
    out = make_mirscore_out(tmp_path, "strucVis", ".ps.pdf")
    for verdict in ("passed", "failed"):
        (out / "RNAplots" / verdict).mkdir(parents=True)
    (out / "RNAplots" / "passed" / f"{P1}_plot.pdf").write_text("worse plot")
    dest = tmp_path / "at.mirna_structures"
    collect_structures(out, recs(), dest)
    assert (dest / "passed" / f"{P1}.ps.pdf").exists()
    assert not (dest / "passed" / f"{P1}_plot.pdf").exists()


def test_no_plot_source_is_a_warning_not_a_crash(tmp_path):
    dest = tmp_path / "at.mirna_structures"
    assert collect_structures(tmp_path / "mirscore", recs(), dest) == 0


def test_a_plot_from_a_previous_run_does_not_survive_this_one(tmp_path):
    """A locus that passed then failed would answer twice, once per directory."""
    dest = tmp_path / "at.mirna_structures"
    (dest / "passed").mkdir(parents=True)
    (dest / "passed" / f"{P2}.ps.pdf").write_text("stale verdict")

    out = make_mirscore_out(tmp_path, "strucVis", ".ps.pdf")
    collect_structures(out, recs(), dest)
    assert not (dest / "passed" / f"{P2}.ps.pdf").exists()
    assert (dest / "failed" / f"{P2}.ps.pdf").exists()


def test_a_scored_prediction_without_a_plot_is_reported(tmp_path, capsys):
    """The source exists and copies nothing for P2: silence would hide it."""
    out = make_mirscore_out(tmp_path, "strucVis", ".ps.pdf")
    (out / "strucVis" / "failed" / f"{P2}.ps.pdf").unlink()
    assert collect_structures(out, recs(), tmp_path / "at.mirna_structures") == 1
    err = capsys.readouterr().err
    assert "WARN" in err
    assert "1" in err and "2" in err


def test_a_never_scored_prediction_has_no_plot(tmp_path):
    """P3 never reached miRScore, so there is nothing to copy for it."""
    out = make_mirscore_out(tmp_path, "strucVis", ".ps.pdf")
    dest = tmp_path / "at.mirna_structures"
    collect_structures(out, recs(), dest)
    assert not list(dest.rglob(f"{P3}*"))
