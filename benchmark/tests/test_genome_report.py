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
import annotation
import report_cli
import report_genome
import report_html
import report_text
import srnaprofiler
import test_html_report_cli as base

GFF = ("##gff-version 3\n"
       "c1\tsrc\tgene\t60\t400\t.\t+\t.\tID=gene:G1;Name=G1\n"
       "c1\tsrc\tmRNA\t60\t400\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
       "c1\tsrc\texon\t60\t400\t.\t+\t.\tID=e1;Parent=transcript:T1\n")


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
    assert len(gm["g"]) == 1 and gm["g"][0][0] == "G1"
    assert gm["g"][0][5] == [["E", 59, 400]]
    assert "cv" in gm and gm["cv"]["s"] == gm["w"][0]
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
              "g": [["Gname", "+", ws + 100, ws + 3000, 0,
                     [["E", ws + 100, ws + 400], ["C", ws + 200, ws + 400]]] for _ in range(2)],
              "t": [["LTR/Gypsy", "Gypsy-1", "+", ws + 1000, ws + 1500, 1] for _ in range(3)],
              "cv": cv}
        total += len(report_html.gm_json(gm))
    per = total / 3000
    print(f"gm bytes per locus: {per:.0f}")
    assert per < 12000


def _tsv(*rows):
    return "\n".join("\t".join(map(str, r)) for r in rows) + "\n"


def _locus(*specs):
    from types import SimpleNamespace as NS
    return [NS(name="L", chrom="c1", start=300, end=350)]


def _windows(tmp_path, text, name):
    p = tmp_path / name
    p.write_text(text)
    layers = annotation.load_layers(["gene:" + str(p)], {"c1"})
    return report_genome.build_windows(_locus(), layers, 2000, {"c1": 5000})["L"]


def test_e2e_ensembl_gff3(tmp_path):
    t = _tsv(["c1", "s", "gene", 101, 1000, ".", "+", ".", "ID=gene:G1;Name=Gname"],
             ["c1", "s", "mRNA", 101, 1000, ".", "+", ".", "ID=transcript:T1;Parent=gene:G1"],
             ["c1", "s", "exon", 101, 400, ".", "+", ".", "Parent=transcript:T1;ID=e1"],
             ["c1", "s", "exon", 701, 1000, ".", "+", ".", "Parent=transcript:T1;ID=e2"],
             ["c1", "s", "CDS", 201, 400, ".", "+", "0", "Parent=transcript:T1;ID=cds1"])
    g = _windows(tmp_path, t, "e.gff3")["g"]
    assert len(g) == 1 and g[0][0] == "Gname" and g[0][1] == "+" and g[0][2:4] == [100, 1000]
    assert g[0][5] == [["E", 100, 400], ["C", 200, 400], ["E", 700, 1000]]


def test_e2e_ncbi_gff3(tmp_path):
    t = _tsv(["c1", "s", "gene", 101, 1000, ".", "-", ".", "ID=gene-LOC1;Name=LOC1;gene=LOC1"],
             ["c1", "s", "mRNA", 101, 1000, ".", "-", ".", "ID=rna-XM_1;Parent=gene-LOC1;gene=LOC1"],
             ["c1", "s", "exon", 101, 400, ".", "-", ".", "ID=exon-XM_1-1;Parent=rna-XM_1"],
             ["c1", "s", "exon", 701, 1000, ".", "-", ".", "ID=exon-XM_1-2;Parent=rna-XM_1"])
    g = _windows(tmp_path, t, "n.gff3")["g"]
    assert len(g) == 1 and g[0][0] == "LOC1" and g[0][1] == "-"
    assert g[0][5] == [["E", 100, 400], ["E", 700, 1000]]


def test_e2e_gtf(tmp_path):
    a = 'gene_id "G1"; gene_name "Gname"; transcript_id "T1";'
    t = _tsv(["c1", "s", "gene", 101, 1000, ".", "+", ".", 'gene_id "G1"; gene_name "Gname";'],
             ["c1", "s", "transcript", 101, 1000, ".", "+", ".", a],
             ["c1", "s", "exon", 101, 400, ".", "+", ".", a],
             ["c1", "s", "exon", 701, 1000, ".", "+", ".", a])
    g = _windows(tmp_path, t, "t.gtf")["g"]
    assert len(g) == 1 and g[0][0] == "Gname"
    assert g[0][5] == [["E", 100, 400], ["E", 700, 1000]]


def test_gm_drops_mo_when_nothing_capped_and_keeps_it_otherwise():
    g = {"w": [0, 9], "L": [1, 2], "g": [], "t": [], "mo": {"g": 0, "t": 0}, "ly": [],
         "cv": {"s": 0}}
    assert report_html._gm(g) == {"w": [0, 9], "L": [1, 2], "cv": {"s": 0}}
    g["mo"] = {"g": 0, "t": 7}
    g["ly"] = ["gene"]
    assert report_html._gm(g)["mo"] == {"t": 7} and report_html._gm(g)["ly"] == ["gene"]


def test_gm_payload_has_no_empty_mo(tmp_path):
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)],
                                        flank=2000))
    p = payload(tmp_path)
    for d in p["details"].values():
        assert "mo" not in d["gm"]
        assert d["gm"]["ly"] == ["gene"]


def test_contract_documents_row_layouts():
    doc = report_html.__doc__
    for frag in ("mo {track: n rows not shown}", "[label, strand, s, e, hl, blocks]",
                 "[label, family, strand, s, e, hl]", "[label, ftype, strand, s, e, hl]",
                 "blocks [t, s, e]"):
        assert frag in doc


def test_flank_zero(tmp_path):
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)],
                                        flank=0))
    p = payload(tmp_path)
    idx = [i for i, r in enumerate(p["rows"]) if r["n"] == "cl1"][0]
    gm = p["details"][str(idx)]["gm"]
    assert gm["w"] == gm["L"] == [100, 200] and p["meta"]["flank"] == 0


def test_contig_absent_from_annotation(tmp_path):
    report_cli.run_report(base.args_for(tmp_path, annotation=["gene:" + gff(tmp_path)],
                                        flank=2000))
    p = payload(tmp_path)
    idx = [i for i, r in enumerate(p["rows"]) if r["n"] == "cl2"][0]   # c2: no features
    gm = p["details"][str(idx)]["gm"]
    assert not any(k in gm for k in ("g", "t", "r", "p", "mo"))
    assert gm["ly"] == ["gene"] and "cv" in gm


def test_plain_report_exact_keys(tmp_path):
    report_cli.run_report(base.args_for(tmp_path))
    p = payload(tmp_path)
    assert set(p["meta"]) == {"title", "floor", "n_total", "kingdom", "mode", "version",
                              "github", "tagline", "install", "defining", "conf_min",
                              "conf_k", "dot_names", "dot_legend", "kt", "n_hidden",
                              "capped", "n_detail", "cl", "reads_shown", "summary"}
    for d in p["details"].values():
        assert set(d) <= {"v", "strip", "ctx", "seq", "rpm", "sig", "an", "ck", "mi"}
        assert not set(d) & {"gm", "flank", "ly", "mo"}


# ui.gm.* keys that no code reads yet. The genome-view JS tasks must REMOVE each key they use
# from this set (the test below fails otherwise). Task 5+ : when app.js reads a key, delete it
# here; a key built dynamically (e.g. "ui.gm.track." + kind) needs the prefix handled in
# test_html_report_html.test_strings_have_no_orphans instead.
# (the figure, Task 5, now reads: title, track.cv/ends/g/t/r/p, empty, nolayer, more, tip.*,
# locus, range, legend.te, legend.size; the zoom/pan keys and track.orf are still unread)
GM_UNUSED = {f"ui.gm.{k}" for k in GM_KEYS} - {f"ui.gm.{k}" for k in (
    "title", "track.cv", "track.ends", "track.g", "track.t", "track.r", "track.p", "empty", "nolayer",
    "more", "tip.type", "tip.strand", "tip.pos", "tip.overlap", "tip.nostrand", "locus", "range",
    "legend.te", "legend.size")}


def test_gm_unused_allowlist_shrinks_with_the_js():
    app = open(os.path.join(os.path.dirname(__file__), "..", "..", "report_assets",
                            "app.js"), encoding="utf-8").read()
    used = set(re.findall(r'"(ui\.gm\.[A-Za-z0-9_.]+)"', app))
    stale = used & GM_UNUSED
    assert not stale, f"app.js uses {sorted(stale)}: remove them from GM_UNUSED"
    assert GM_UNUSED <= set(report_text.STRINGS)
