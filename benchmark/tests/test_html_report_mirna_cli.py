"""--mirna on the report subcommand: auto-detect, explicit path, matching, warnings."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.dirname(__file__))
import bamutil
import report_cli
import report_mirna as rm
from test_html_report_cli import COLS, args_for
from test_html_report_html import data_of
from test_html_report_mirna import (END, HI, LO, M_AB, START, plus_reads, precursor, row_dict,
                                    write_tsv)


def setup(tmp_path, cls="miRNA", mirna_name="run.mirna.tsv", **kw):
    a = args_for(tmp_path, **kw)            # writes its own inputs, replaced below
    clusters = tmp_path / "run.clusters.tsv"
    rows = [["cl1", "c1", START, END, 500, 50, cls, 1.0, 22, 0.8, 0.5, 0.8, 0.5, 3, 0.5, 0.1, 0.9,
             "False", "UGA", ""],
            ["cl2", "c1", 3000, 3100, 400, 40, "piRNA", 0.9, 26, 0.8, 0.5, 0.8, 0.5, 3, 0.5, 0.1, 0.9,
             "False", "UGA", ""]]
    clusters.write_text("\n".join("\t".join(map(str, r)) for r in [COLS] + rows) + "\n")
    bamutil.make_seq_bam(a.bam, plus_reads(precursor()) + [("c1", 3000, "A" * 26, "+", 1)])
    for old in tmp_path.glob("*.mirna.tsv"):
        old.unlink()
    if mirna_name:
        write_tsv(tmp_path / mirna_name, [row_dict()])
    return a


def details(tmp_path):
    out = tmp_path / "out" / "x" / "report.report.html"
    d = data_of(out.read_text())
    return {r["n"]: d["details"][str(i)] for i, r in enumerate(d["rows"]) if r["d"]}


def test_default_mirna_file_beside_clusters(tmp_path, monkeypatch):
    monkeypatch.setattr(rm.shutil, "which", lambda n: None)
    monkeypatch.setattr(rm, "_import_rna", lambda: None)
    a = setup(tmp_path)
    assert report_cli._mirna_path(a) == str(tmp_path / "run.mirna.tsv")
    report_cli.run_report(a)
    det = details(tmp_path)
    mi = det["cl1"]["mi"]
    assert mi["ns"] == "m" and "db" not in mi and mi["m"] == list(M_AB) and mi["lo"] == LO
    assert mi["rd"] and len(mi["cov"]) == HI - LO and mi["anti"] == 1.0
    assert "mi" not in det["cl2"]                                # the piRNA cluster has none


def test_explicit_path_and_missing_file(tmp_path):
    a = setup(tmp_path, mirna_name=None)
    assert report_cli._mirna_path(a) is None
    report_cli.run_report(a)                                     # absent auto file: no figure, no error
    assert "mi" not in details(tmp_path)["cl1"]
    other = tmp_path / "elsewhere.tsv"
    write_tsv(other, [row_dict()])
    a = setup(tmp_path, mirna_name=None, mirna=str(other))
    report_cli.run_report(a)
    assert "mi" in details(tmp_path)["cl1"]
    a = setup(tmp_path, mirna_name=None, mirna=str(tmp_path / "nope.tsv"))
    with pytest.raises(SystemExit, match="miRNA TSV not found"):
        report_cli.run_report(a)


def test_only_miRNA_clusters_are_matched(tmp_path):
    report_cli.run_report(setup(tmp_path, cls="siRNA"))
    assert "mi" not in details(tmp_path)["cl1"]


def test_unmatched_row_is_ignored(tmp_path):
    a = setup(tmp_path, mirna_name=None)
    write_tsv(tmp_path / "run.mirna.tsv", [row_dict(chrom="c2")])
    report_cli.run_report(a)
    assert "mi" not in details(tmp_path)["cl1"]


def test_bad_table_is_an_error(tmp_path):
    a = setup(tmp_path)
    (tmp_path / "run.mirna.tsv").write_text("x\ty\n1\t2\n")
    with pytest.raises(SystemExit, match="lacks the column"):
        report_cli.run_report(a)


def test_warns_when_rnafold_is_missing(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(rm.shutil, "which", lambda n: None)
    monkeypatch.setattr(rm, "_import_rna", lambda: None)
    report_cli.run_report(setup(tmp_path))
    err = capsys.readouterr().err
    assert "RNAfold not found" in err and "1 of 1 miRNA clusters matched" in err


def test_flag_is_on_the_subcommand():
    import subprocess
    out = subprocess.run([sys.executable, os.path.join(os.path.dirname(__file__), "..", "..",
                          "srnaprofiler.py"), "report", "--help"], capture_output=True, text=True)
    assert "--mirna" in out.stdout
