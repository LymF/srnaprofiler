"""Task 4: report --annotation / --flank, genome-view payload and contract."""
import gzip
import json
import os
import re
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.dirname(__file__))
import report_cli
import report_html
import report_text
import srnaprofiler
import test_html_report_cli as base

GFF = ("##gff-version 3\n"
       "c1\tsrc\tgene\t60\t400\t.\t+\tID=gene:G1;Name=G1\n"
       "c1\tsrc\tmRNA\t60\t400\t.\t+\tID=transcript:T1;Parent=gene:G1\n"
       "c1\tsrc\texon\t60\t400\t.\t+\tID=e1;Parent=transcript:T1\n")


def payload(tmp_path):
    out = tmp_path / "out" / "x" / "report.report.html"
    m = re.search(r'<script type="application/json" id="report-data">(.*?)</script>',
                  out.read_text(), re.S)
    return json.loads(m.group(1))


def gff(tmp_path):
    p = tmp_path / "g.gff3"
    p.write_text(GFF)
    return str(p)


def parse(*extra):
    return srnaprofiler.build_parser().parse_args(
        ["report", "--clusters", "a", "--bam", "b", "-o", "o", *extra])


def test_parser_flags():
    a = parse()
    assert a.annotation is None and a.flank == 2000
    a = parse("--annotation", "gene:x.gff", "--annotation", "te:y.out", "--flank", "500")
    assert a.annotation == ["gene:x.gff", "te:y.out"] and a.flank == 500


def test_negative_flank_errors(tmp_path):
    a = base.args_for(tmp_path, annotation=None, flank=-1)
    with pytest.raises(SystemExit):
        report_cli.run_report(a)


def test_plain_report_has_no_gm(tmp_path):
    report_cli.run_report(base.args_for(tmp_path))
    p = payload(tmp_path)
    assert all("gm" not in d for d in p["details"].values())
    assert "gm" not in p["meta"] and "flank" not in p["meta"]
    assert '"gm"' not in json.dumps(p)


def test_plain_payload_identical_with_flags_absent_or_none(tmp_path):
    report_cli.run_report(base.args_for(tmp_path))
    one = payload(tmp_path)
    report_cli.run_report(base.args_for(tmp_path, annotation=None, flank=2000))
    assert payload(tmp_path) == one


def test_gm_in_detail(tmp_path):
    a = base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=2000)
    report_cli.run_report(a)
    p = payload(tmp_path)
    assert p["meta"]["gm"] == 1 and p["meta"]["flank"] == 2000
    idx = [i for i, r in enumerate(p["rows"]) if r["n"] == "cl1"][0]
    gm = p["details"][str(idx)]["gm"]
    assert gm["w"] == [max(0, 100 - 2000), min(5000, 200 + 2000)]
    assert gm["L"] == [100, 200]
    assert gm["g"] and "cv" in gm and gm["cv"]["s"] == gm["w"][0]
    assert gm["cv"]["nb"] <= 200


def test_flank_changes_window(tmp_path):
    a = base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)], flank=500)
    report_cli.run_report(a)
    p = payload(tmp_path)
    idx = [i for i, r in enumerate(p["rows"]) if r["n"] == "cl1"][0]
    assert p["details"][str(idx)]["gm"]["w"] == [0, 700]
    assert p["meta"]["flank"] == 500


def _spy_bam(monkeypatch):
    calls = []
    monkeypatch.setattr(report_cli.report_signals, "signals_for",
                        lambda *a, **k: calls.append(1) or {})
    monkeypatch.setattr(report_cli.report_signals, "library_lengths",
                        lambda *a, **k: calls.append(1) or [])
    return calls


def test_bad_label_exits_before_bam_work(tmp_path, monkeypatch):
    calls = _spy_bam(monkeypatch)
    a = base.args_for(tmp_path, annotation=["bogus:" + gff(tmp_path)], flank=2000)
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert e.value.code not in (0, None) and not calls


def test_missing_file_exits_before_bam_work(tmp_path, monkeypatch):
    calls = _spy_bam(monkeypatch)
    a = base.args_for(tmp_path, annotation=["gene:" + str(tmp_path / "no.gff")], flank=2000)
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert e.value.code not in (0, None) and not calls


def test_fake_gz_exits_before_bam_work(tmp_path, monkeypatch):
    calls = _spy_bam(monkeypatch)
    f = tmp_path / "fake.gff3.gz"
    f.write_bytes(b"this is not gzip")
    a = base.args_for(tmp_path, annotation=["gene:" + str(f)], flank=2000)
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert e.value.code not in (0, None) and not calls


def test_truncated_gz_exits(tmp_path, monkeypatch):
    calls = _spy_bam(monkeypatch)
    f = tmp_path / "t.gff3.gz"
    f.write_bytes(gzip.compress((GFF * 200).encode())[:-20])
    a = base.args_for(tmp_path, annotation=["gene:" + str(f)], flank=2000)
    with pytest.raises(SystemExit) as e:
        report_cli.run_report(a)
    assert e.value.code not in (0, None) and not calls


def test_contract_mentions_gm():
    doc = report_html.__doc__
    assert "gm" in doc and "meta.gm" in doc and "flank" in doc


GM_KEYS = ["title", "track.cv", "track.ends", "track.g", "track.t", "track.r", "track.p",
           "track.orf", "empty", "nolayer", "more", "tip.type", "tip.strand", "tip.pos",
           "tip.overlap", "tip.nostrand", "zoom.in", "zoom.out", "zoom.reset", "pan.left",
           "pan.right", "locus", "range", "legend.te", "legend.size"]


def test_gm_strings():
    for k in GM_KEYS:
        assert report_text.STRINGS[f"ui.gm.{k}"]
    S = report_text.STRINGS
    assert "{n}" in S["ui.gm.more"] and "{kind}" in S["ui.gm.nolayer"]
    assert all(x in S["ui.gm.range"] for x in ("{chrom}", "{start}", "{end}", "{n}"))


def test_payload_size_per_locus():
    """Realistic synthetic gm: 3000 loci, a few features each; logs bytes per locus."""
    import random
    rnd = random.Random(1)
    from report_signals import _window_build
    from collections import Counter
    total = 0
    for i in range(3000):
        ws = rnd.randrange(0, 10**7)
        keys = Counter()
        for _ in range(300):
            keys[("+" if rnd.random() < .5 else "-", rnd.randrange(ws, ws + 4380),
                  rnd.choice((21, 22, 24, 26)), rnd.choice((1, 1, 2, 4)), "A")] += 1
        cv = _window_build(keys, ws, ws + 4400, 200)
        gm = {"w": [ws, ws + 4400], "L": [ws + 2000, ws + 2100],
              "g": [["Gname", "gene", "+", ws + 100, ws + 3000, [[ws + 100, ws + 400, "exon"]], 0]
                    for _ in range(2)],
              "t": [["LTR/Gypsy", "Gypsy-1", "+", ws + 1000, ws + 1500, 1] for _ in range(3)],
              "cv": cv}
        total += len(report_html.gm_json(gm))
    per = total / 3000
    print(f"gm bytes per locus: {per:.0f}")
    assert per < 12000
