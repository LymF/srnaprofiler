"""Task 8: real-data fixes (A1, A2) and the visual polish helpers (B1-B10)."""
import json
import logging
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import test_genome_figure_js as t5  # noqa: E402
from annotation import Feature, IntervalIndex, Layer, load_layers, summarize_locus  # noqa: E402
from test_genome_figure_js import data_with, gm_payload  # noqa: E402
from test_html_report_js import APP, run_js  # noqa: E402

# ---------------- A1, A2 (plain Python) ----------------
DM6 = (
    "##gff-version 3\n"
    "2L\tFlyBase\tgene\t101\t500\t.\t+\t.\tID=FBgn1;Name=HA8\n"
    "2L\tFlyBase\tmRNA\t101\t500\t.\t+\t.\tID=FBtr1;Parent=FBgn1;Name=HA8-RA\n"
    "2L\tFlyBase\texon\t101\t500\t.\t+\t.\tParent=FBtr1\n"
    "2L\tFlyBase\ttransposable_element_gene\t601\t900\t.\t+\t.\tID=FBgn2;Name=invader1{}2621\n"
    "2L\tFlyBase\ttransposable_element\t601\t900\t.\t+\t.\tID=FBti2;Name=invader1{}2621\n"
    "2L\tFlyBase\tTransposable_Element_Gene\t901\t950\t.\t+\t.\tID=FBgn3;Name=X\n")


def _layer(tmp_path, text):
    p = tmp_path / "dm6.gff3"
    p.write_text(text)
    return load_layers([f"gene:{p}"], {"2L"})[0]


def test_a1_te_insertions_never_reach_the_gene_layer(tmp_path, caplog):
    caplog.set_level(logging.DEBUG)
    L = _layer(tmp_path, DM6)
    kept = sorted(f.ftype for f in L.index["2L"].overlaps(0, 2000))
    assert kept == ["exon", "gene", "mRNA"]
    assert L.n_features == 3
    # a TE insertion is intergenic for the gene layer, not a gene_intron
    assert summarize_locus("2L", 610, 700, "+", [L], 0.5)["contexts"] == {"intergenic"}


def test_a1_ignored_count_includes_the_te_rows(tmp_path, capsys):
    _layer(tmp_path, DM6)
    cap = capsys.readouterr()
    assert "ignored 3 features" in (cap.out + cap.err)


def test_a1_class_function_and_other_layers_unaffected():
    from annotation import _gene_class
    for t in ("transposable_element_gene", "Transposable_Element", "transposable_element"):
        assert _gene_class(t) is None
    assert _gene_class("gene") == "span" and _gene_class("ncRNA_gene") == "span"
    L = Layer("te", "x", {"c": IntervalIndex([Feature(0, 100, "+", "transposable_element", "LTR/Gypsy", "t1")])}, 1)
    assert summarize_locus("c", 10, 30, "+", [L], 0.5)["text"] == "te:LTR/Gypsy:sense:1.00:t1"


def test_a2_gene_entry_uses_the_label_when_present(tmp_path):
    L = _layer(tmp_path, DM6)
    r = summarize_locus("2L", 150, 300, "+", [L], 0.5)
    assert r["text"].startswith("gene:HA8:sense:") or r["text"].startswith("gene:HA8-RA:sense:") \
        or r["text"].startswith("gene:FBtr1:sense:") or r["text"].startswith("gene:exon:sense:")
    G = Layer("gene", "x", {"c": IntervalIndex([Feature(0, 1000, "+", "gene", "gene", "g1", label="HA8")])}, 1)
    assert summarize_locus("c", 10, 30, "+", [G], 0.5)["text"] == "gene:HA8:sense:1.00:g1"
    G2 = Layer("gene", "x", {"c": IntervalIndex([Feature(0, 1000, "+", "gene", "gene", "g1")])}, 1)
    assert summarize_locus("c", 10, 30, "+", [G2], 0.5)["text"] == "gene:gene:sense:1.00:g1"
    T = Layer("te", "x", {"c": IntervalIndex([Feature(0, 100, "+", "feature", "LTR/Gypsy", "t1", label="TE_9")])}, 1)
    assert summarize_locus("c", 10, 30, "+", [T], 0.5)["text"] == "te:LTR/Gypsy:sense:1.00:t1"


def test_a2_exon_driven_entry_names_the_gene(tmp_path):
    # the exon row has no Name of its own: the entry still carries the gene's name
    L = _layer(tmp_path, DM6)
    assert summarize_locus("2L", 150, 300, "+", [L], 0.5)["text"].split(":")[1] == "HA8"


# ---------------- B (stub DOM) ----------------
pytestmark_js = pytest.mark.skipif(shutil.which("osascript") is None, reason="needs osascript (JavaScriptCore)")

HOOK = ("window.__gm = {gmLanes: gmLanes, gmClip: gmClip, gmCoverage: gmCoverage, gmTeClass: gmTeClass,"
        " gmPlace: gmPlace, gmTrackH: gmTrackH, tipPos: tipPos, gmSigned: gmSigned};\n  buildTabs();\n  renderAbout();")
APPX = APP.replace("  buildTabs();\n  renderAbout();", HOOK, 1)
assert APPX != APP

O = [[1500, 1800, "+", 99], [3000, 3600, "-", 199]]


def run(tmp_path, d, body="return res;", width=900):
    return run_js(tmp_path, json.dumps(d), t5.OPEN + body, app=APPX, width=width)


def orf_data(tmp_path, o=O, **over):
    gm = gm_payload(**over)
    if o is not None:
        gm["o"] = o
    d = data_with(tmp_path, gm)
    d["meta"]["orf_min_aa"] = 100
    return d


@pytestmark_js
def test_b6_button_is_named_zoom_to_locus_and_band_keeps_locus(tmp_path):
    d = data_with(tmp_path, gm_payload())
    S = d["strings"]
    assert S["ui.gm.locus"] == "Zoom to locus" and S["ui.gm.band"] == "Locus"
    body = """
    res.btns = H.all(gmEl, "button").map(H.text);
    res.texts = H.all(gmEl, "text").map(H.text);
    return res;
    """
    r = run(tmp_path, d, body)
    assert "Zoom to locus" in r["btns"] and "Locus" not in r["btns"]
    assert "Locus" in r["texts"]


@pytestmark_js
def test_b6_controls_are_grouped_and_disabled_is_readable(tmp_path):
    body = """
    var g = H.byClass(gmEl, "rpt-gm-grp", "span");
    return g.map(function (x) { return H.all(x, "button").length; });
    """
    assert run(tmp_path, data_with(tmp_path, gm_payload()), body) == [2, 2, 2]
    css = open(os.path.join(os.path.dirname(__file__), "..", "..", "report_assets", "report.css"), encoding="utf-8").read()
    import re
    m = re.search(r"\.rpt-gm-ctl \.rpt-btn\[disabled\]\s*\{([^}]*)\}", css)
    assert m and float(re.search(r"opacity:\s*([\d.]+)", m.group(1)).group(1)) >= 0.55


@pytestmark_js
def test_b1_orf_label_short_caption_and_height(tmp_path):
    d = orf_data(tmp_path)
    S = d["strings"]
    assert S["ui.gm.track.orf"] == "Predicted ORFs" and S["ui.gm.orf.note"] == "experimental"
    body = """
    var cap = H.all(gmEl, "text").filter(function (t) { return /rpt-gm-orfnote/.test(H.cls(t)); }).map(H.text);
    var tl = H.all(gmEl, "text").filter(function (t) { return H.cls(t) === "rpt-gm-tl"; }).map(H.text);
    res.cap = cap; res.tl = tl;
    res.h = {empty: window.__gm.gmTrackH("o", "empty", 0, 12, 0), geneEmpty: window.__gm.gmTrackH("g", "empty", 0, 26, 0),
             teEmpty: window.__gm.gmTrackH("t", "empty", 0, 12, 0), orfRows: window.__gm.gmTrackH("o", "", 2, 12, 0),
             geneRows: window.__gm.gmTrackH("g", "", 3, 26, 0), more: window.__gm.gmTrackH("t", "", 1, 12, 5)};
    return res;
    """
    r = run(tmp_path, d, body)
    assert r["cap"] == ["experimental"] and "Predicted ORFs" in r["tl"]
    h = r["h"]
    assert h["empty"] == h["geneEmpty"] == h["teEmpty"]
    assert h["orfRows"] >= h["empty"] and h["orfRows"] <= 40
    assert h["geneRows"] == 3 * 26 + 12 and h["more"] >= 46


@pytestmark_js
def test_b2_antisense_ticks_use_a_real_minus_and_ends_have_a_unit(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    res.tk = H.all(gmEl, "text").filter(function (t) { return H.cls(t) === "tk-y"; }).map(H.text);
    res.ts = H.all(gmEl, "text").filter(function (t) { return /rpt-gm-ts/.test(H.cls(t)); }).map(H.text);
    res.s = [window.__gm.gmSigned(5, 1), window.__gm.gmSigned(5, -1), window.__gm.gmSigned("0", 1)];
    return res;
    """
    r = run(tmp_path, d, body)
    assert r["s"] == ["5", "−5", "0"]
    neg = [t for t in r["tk"] if t.startswith("−")]
    assert len(neg) >= 2 and all("-" not in t for t in r["tk"])
    assert d["strings"]["ui.gm.ends.unit"] in " ".join(r["ts"])


@pytestmark_js
def test_b4_legend_groups_keep_heading_with_entries_and_has_mark_entries(tmp_path):
    d = data_with(tmp_path, gm_payload(), mi={"n": "x", "lo": 2900, "L": 120, "nowin": 0, "t": [], "pc": []})
    S = d["strings"]
    body = """
    var ul = H.byClass(gmEl, "rpt-srna-legend", "ul")[0];
    res.txt = H.text(ul);
    res.grp = H.byClass(ul, "grp", "li").map(function (li) {
      var kids = li.children || li.childNodes;
      return {first: H.text(kids[0]), n: kids.length, firstEntries: H.all(kids[0], "li").length + H.all(kids[0], "span").length};
    });
    res.heads = H.byClass(ul, "hd").length;
    return res;
    """
    r = run(tmp_path, d, body)
    for k in ("ui.gm.legend.overlap", "ui.gm.legend.band", "ui.gm.legend.cont"):
        assert S[k] in r["txt"], k
    assert S["mi.title"] in r["txt"]
    assert r["grp"] and all(g["n"] >= 2 for g in r["grp"])
    assert any(S["ui.gm.legend.te"] in g["first"] for g in r["grp"])
    # the unit that holds a heading also holds its first entry
    for g in r["grp"]:
        assert len(g["first"].strip()) > len(g["first"].split()[0]) or g["firstEntries"] > 1


@pytestmark_js
def test_b9_figure_has_no_letter_but_keeps_its_hook(tmp_path):
    body = """
    res.letters = H.byClass(gmEl, "rpt-srna-fig__letter").length;
    res.title = H.text(H.all(gmEl, "h3")[0]);
    return res;
    """
    d = data_with(tmp_path, gm_payload())
    r = run(tmp_path, d, body)
    assert r["has"] and r["letters"] == 0 and r["title"].startswith(d["strings"]["ui.gm.title"])


@pytestmark_js
def test_b8_tooltip_position_stays_inside_the_viewport(tmp_path):
    body = """
    var P = window.__gm.tipPos, out = [];
    out.push(P(100, 50, 200, 80, 1000, 600));          // plenty of room: below right of the pointer
    out.push(P(100, 580, 200, 80, 1000, 600));         // near the bottom: flips above
    out.push(P(980, 300, 200, 80, 1000, 600));         // near the right edge: flips left
    out.push(P(100, 300, 200, 700, 1000, 600));        // taller than the viewport: pinned inside
    out.push(P(100, 20, 200, 560, 1000, 600));         // fits neither way: takes the roomier side, never off screen
    return out;
    """
    r = run(tmp_path, data_with(tmp_path, gm_payload()), body)
    assert r[0] == [114, 64] and r[1][1] == 580 - 14 - 80 and r[2][0] == 980 - 14 - 200
    for (x, y), (w, h) in zip(r[1:], [(200, 80), (200, 80), (200, 700), (200, 560)]):
        assert x >= 8 and y >= 8
    assert r[4][1] + 560 <= 600 - 8
    assert r[3][1] == 8


@pytestmark_js
def test_b3_label_placement_rules(tmp_path):
    body = """
    function xa(p) { return 112 + (p - 1000) / 5000 * 772; }
    function row(blocks, s, e, label) { return {kind: "g", label: label, strand: "+", s: s, e: e, hl: false, blocks: blocks}; }
    var long = "LongGeneName12";                    // 14 chars: 112 px
    var wide = row([["E", 1000, 1100], ["E", 5900, 6000]], 1000, 6000, long);       // whole window, wide gap between two exons
    var solid = row([["E", 1000, 6000]], 1000, 6000, long);                        // one exon over everything
    var narrow = row([["E", 1000, 3000], ["E", 3050, 6000]], 1000, 6000, long);     // gap too small for the label
    var near = row([["E", 3000, 3100]], 3000, 3100, "G1");                          // room beside
    var out = window.__gm.gmPlace("g", [wide, solid, narrow, near], xa, 112, 884, 1000, 6000).items.map(function (it) {
      return it.label ? {mode: it.label.mode, x: it.label.x, anchor: it.label.anchor, w: it.label.w} : null;
    });
    return out;
    """
    r = run(tmp_path, data_with(tmp_path, gm_payload()), body)
    assert r[0]["mode"] == "pill" and r[1]["mode"] == "box" and r[2]["mode"] == "box" and r[3]["mode"] == "beside"
    for lb in r:
        half = lb["w"] / 2 if lb["anchor"] == "middle" else 0
        lo = lb["x"] - half if lb["anchor"] != "end" else lb["x"] - lb["w"]
        hi = lo + lb["w"]
        assert lo >= 112 - 1e-6 and hi <= 884 + 1e-6, lb


@pytestmark_js
def test_b7_zoomed_locus_label_stays_inside_the_plot(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    H.click(H.buttons(gmEl, %s)[0]); H.frames();
    var svg = H.all(gmEl, "svg")[0], W = +svg.getAttribute("width");
    var t = H.all(gmEl, "text").filter(function (e) { return H.cls(e) === "rpt-gm-loclabel"; })[0];
    return {x: +t.getAttribute("x"), a: t.getAttribute("text-anchor"), W: W, txt: H.text(t)};
    """ % json.dumps(d["strings"]["ui.gm.locus"])
    r = run(tmp_path, d, body)
    assert r["txt"] == "Locus"
    w = 5 * 8
    lo = r["x"] - (w if r["a"] == "end" else 0)
    assert lo >= 112 and lo + w <= r["W"] - 16 + 1e-6
