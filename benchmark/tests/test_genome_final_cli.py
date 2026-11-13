"""Final-review fixes in report: size guard, NA clusters, flank bounds, ignored flags, no-contig layers."""
import os
import sys
from collections import namedtuple

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import report_cli  # noqa: E402
import test_html_report_cli as base  # noqa: E402
from annotation import Feature, IntervalIndex, Layer  # noqa: E402
from report_genome import build_windows, window_features  # noqa: E402
from test_genome_report import gff, payload  # noqa: E402


def _err(capsys):
    c = capsys.readouterr()
    return c.out + c.err


def test_size_guard_compares_total_and_uses_mib(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(report_cli, "_SIZE_WARN_BYTES", 10)
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=2000))
    t = _err(capsys)
    assert "[WARN]" in t and "MiB" in t and "--max-detail" in t and "--flank" not in t.split("[WARN]")[-1]
    assert "adds" in t and "MiB" in t


def test_size_guard_silent_below_threshold(tmp_path, capsys):
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=2000))
    assert "[WARN]" not in _err(capsys)


def test_size_guard_applies_without_annotation(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(report_cli, "_SIZE_WARN_BYTES", 10)
    report_cli.run_report(base.args_for(tmp_path))
    assert "MiB" in _err(capsys)


def test_na_start_end_skipped_with_one_warn(capsys):
    C = namedtuple("C", "name chrom start end")
    out = build_windows([C("a", "c", None, 20), C("b", "c", 10, None), C("ok", "c", 10, 20)],
                        [], 5, {"c": 100})
    assert set(out) == {"ok"}
    t = _err(capsys)
    assert t.count("[WARN]") == 1 and "2 clusters" in t


def test_flank_upper_bound(tmp_path):
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=100001))
    assert "--flank" in str(e.value) and "100000" in str(e.value)
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=100000))


def test_ignored_flags_warn_once_each(tmp_path, capsys):
    report_cli.run_report(base.args_for(tmp_path, flank=500, genome="x.fa", orf_min_aa=50))
    t = _err(capsys)
    for flag in ("--flank", "--genome", "--orf-min-aa"):
        assert t.count(f"{flag} is ignored") == 1, flag
    assert t.count("[WARN]") == 3


def test_flags_given_with_their_partner_do_not_warn(tmp_path, capsys):
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=500))
    assert "ignored" not in _err(capsys)


def test_defaults_do_not_warn(tmp_path, capsys):
    report_cli.run_report(base.args_for(tmp_path, flank=2000, genome=None, orf_min_aa=100))
    assert "[WARN]" not in _err(capsys)


def test_layer_without_matching_contig_is_nm_not_ly():
    empty = Layer("te", "x", {}, 0)
    gene = Layer("gene", "x", {"c": IntervalIndex([Feature(0, 5, "+", "gene", "g")])}, 1)
    out = window_features("c", 0, 100, (10, 20), [empty, gene])
    assert out["ly"] == ["gene"] and out["nm"] == ["te"]
    out = window_features("c", 0, 100, (10, 20), [gene])
    assert "nm" not in out


def test_nm_reaches_the_payload_and_is_dropped_when_empty(tmp_path):
    te = tmp_path / "t.bed"
    te.write_text("chrZ\t1\t5\tx\n")
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path), "te:" + str(te)], flank=2000))
    gms = [d["gm"] for d in payload(tmp_path)["details"].values() if "gm" in d]
    assert gms and all(g["nm"] == ["te"] and g["ly"] == ["gene"] for g in gms)
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=2000))
    assert all("nm" not in d["gm"] for d in payload(tmp_path)["details"].values() if "gm" in d)


# ---- I-4: classify and report annotation consistency ----

def _run_json(tmp_path, ann):
    import json
    d = {"min_reads": 1}
    if ann is not None:
        d["annotation"] = ann
    (tmp_path / "run.run.json").write_text(json.dumps(d))


def test_consistency_warns_on_difference(tmp_path, capsys):
    _run_json(tmp_path, ["gene:/other/dir/genes.gff3"])
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=2000))
    t = _err(capsys)
    assert t.count("differs from the classify run") == 1 and "genes.gff3" in t and "g.gff3" in t


def test_consistency_warns_when_one_side_is_empty(tmp_path, capsys):
    _run_json(tmp_path, None)
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=2000))
    assert "classify used none" in _err(capsys)
    _run_json(tmp_path, ["gene:x.gff3"])
    report_cli.run_report(base.args_for(tmp_path))
    assert "report uses none" in _err(capsys)


def test_consistency_silent_when_same_set_or_no_run_json(tmp_path, capsys):
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=2000))
    assert "differs" not in _err(capsys)
    _run_json(tmp_path, ["gene:/elsewhere/g.gff3"])    # same label and basename
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=2000))
    assert "differs" not in _err(capsys)
    (tmp_path / "run.run.json").write_text("not json")
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=2000))
    assert "differs" not in _err(capsys)
