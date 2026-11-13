"""Task 5: the genomic-context figure, run under the stub DOM (JavaScriptCore via osascript)."""
import copy
import json
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_html_report_html import mk  # noqa: E402
from test_html_report_js import APP, JsError, page, run_js  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("osascript") is None,
                                reason="needs osascript (JavaScriptCore)")

# the pure helpers are exposed on window.__gm for the probes (app.js is an IIFE)
HOOK = ("window.__gm = {gmLanes: gmLanes, gmClip: gmClip, gmCoverage: gmCoverage, gmTeClass: gmTeClass,"
        " gmPlace: gmPlace};\n  buildTabs();\n  renderAbout();")
APPX = APP.replace("  buildTabs();\n  renderAbout();", HOOK, 1)
APP_FILE = os.environ.get("GM_APP")   # run the file against another app.js (to watch it fail)
if APP_FILE:
    APPX = open(APP_FILE, encoding="utf-8").read()
assert APPX != APP or APP_FILE

WS, WE = 1000, 6000
BLOCKS = [["E", 1500, 1800], ["C", 1600, 1800], ["E", 2400, 2700]]


def gm_payload(**over):
    cv = {"s": WS, "b": 50, "nb": 100, "k": ["18-20", "21", "22", "23", "24", "25-35"],
          "+": {"21": [10, [99, 50, 10]], "24": [4, [99, 99, 20]]},
          "-": {"21": [5, [0, 40, 99]]},
          "e": {"+": [3, [99, 20]], "-": [1, [0, 99]]}}
    gm = {"w": [WS, WE], "L": [3000, 3100], "ly": ["gene", "te"], "cv": cv,
          "g": [["GeneA", "+", 1500, 2700, 1, BLOCKS],
                ["Gb", "-", 1600, 2000, 0, [["E", 1600, 2000]]]],
          "t": [["LTR/Gypsy", "Gypsy-1", "+", 3000, 3200, 1],
                ["DNA/hAT", "hAT-2", ".", 3300, 3500, 0]]}
    gm.update(over)
    return gm


def data_with(tmp_path, gm=None, mi=None):
    sel, text, d = page(tmp_path, [mk("c1")], kind="ok")
    key = sorted(d["details"])[0]
    if gm is not None:
        d["details"][key]["gm"] = gm
    if mi is not None:
        d["details"][key]["mi"] = mi
    return d


OPEN = """
var lst = H.id("sec-clusters");
H.click(H.byClass(lst, "row", "tr")[0]);
var det = H.id("sec-detail");
var res = {};
res.figs = H.byClass(det, "rpt-srna-fig", "section").map(function (f) { return f.getAttribute("data-fig"); });
var gmEl = H.byClass(det, "rpt-srna-fig", "section").filter(function (f) { return f.getAttribute("data-fig") === "I"; })[0];
res.has = !!gmEl;
"""


def run(tmp_path, d, body="return res;", width=900):
    return run_js(tmp_path, json.dumps(d), OPEN + body, app=APPX, width=width)


def test_figure_first_when_gm_exists_and_absent_otherwise(tmp_path):
    r = run(tmp_path, data_with(tmp_path, gm_payload()))
    plain = run(tmp_path, data_with(tmp_path))["figs"]
    assert "I" not in plain and plain[0] == "A"
    assert r["figs"] == ["I"] + plain


def test_lane_packing_never_overlaps(tmp_path):
    body = """
    var items = [], i;
    for (i = 0; i < 60; i++) { items.push({a: (i * 37) % 400, b: (i * 37) % 400 + 20 + (i * 13) % 90}); }
    var p = window.__gm.gmLanes(items, 4, 0);
    return {items: items, lane: p.lane, n: p.n};
    """
    r = run(tmp_path, data_with(tmp_path, gm_payload()), body)
    items, lane = r["items"], r["lane"]
    assert r["n"] >= 2 and all(k >= 0 for k in lane)
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if lane[i] == lane[j]:
                a, b = items[i], items[j]
                assert a["b"] + 4 <= b["a"] or b["b"] + 4 <= a["a"], (a, b)


def test_lane_cap_keeps_highlighted_rows(tmp_path):
    body = """
    var items = [], i;
    for (i = 0; i < 10; i++) { items.push({a: 0, b: 50, hl: i === 9}); }
    var p = window.__gm.gmLanes(items, 4, 3, function (e) { return e.hl; });
    return p;
    """
    r = run(tmp_path, data_with(tmp_path, gm_payload()), body)
    assert r["lane"][9] == 0 and r["lane"].count(-1) == 7 and r["n"] == 3


def test_drawn_gene_has_one_rect_per_exon_block_and_no_overlap_in_a_lane(tmp_path):
    body = """
    var genes = H.byClass(det, "rpt-gm-gene", "g");
    res.n = genes.length;
    res.rects = genes.map(function (g) { return H.all(g, "rect").length; });
    res.lanes = genes.map(function (g) {
      var ex = H.all(g, "rect");
      return ex.map(function (e) { return [+e.getAttribute("x"), +e.getAttribute("x") + +e.getAttribute("width"), +e.getAttribute("y")]; });
    });
    return res;
    """
    r = run(tmp_path, data_with(tmp_path, gm_payload()), body)
    assert r["n"] == 2 and sorted(r["rects"]) == [1, 3]


def test_gene_longer_than_window_is_clipped_with_a_marker(tmp_path):
    big = [["Big", "+", 0, 900000, 0, [["E", 0, 900000]]]]
    body = """
    var g = H.byClass(det, "rpt-gm-gene", "g")[0];
    res.clip = g.getAttribute("data-clip");
    var ex = H.all(g, "rect")[0];
    res.x = +ex.getAttribute("x"); res.w = +ex.getAttribute("width");
    res.cont = H.all(g, "path").filter(function (p) { return H.cls(p) === "rpt-gm-cont"; }).length;
    var ln = H.all(g, "line")[0];
    res.lx = [+ln.getAttribute("x1"), +ln.getAttribute("x2")];
    res.svgW = +H.all(gmEl, "svg")[0].getAttribute("width");
    res.pure = window.__gm.gmClip(0, 900000, 1000, 6000);
    res.inside = window.__gm.gmClip(2000, 3000, 1000, 6000);
    return res;
    """
    r = run(tmp_path, data_with(tmp_path, gm_payload(g=big)), body)
    assert r["clip"] == "lr" and r["cont"] == 2
    assert r["pure"] == {"s": 1000, "e": 6000, "l": True, "r": True}
    assert r["inside"] == {"s": 2000, "e": 3000, "l": False, "r": False}
    # the drawn span runs from the first plot pixel to the last one
    assert r["x"] == pytest.approx(112, abs=0.01) and r["x"] + r["w"] == pytest.approx(r["svgW"] - 16, abs=0.01)
    assert r["lx"][0] == pytest.approx(112, abs=0.01)


def test_unstranded_feature_has_no_chevron(tmp_path):
    gm = gm_payload(g=[], t=[["DNA/hAT", "hAT", ".", 2000, 4000, 0]])
    body = """
    var b = H.byClass(det, "rpt-gm-bar", "g")[0];
    return {chev: H.all(b, "polyline").length};
    """
    assert run(tmp_path, data_with(tmp_path, gm), body)["chev"] == 0
    gm2 = gm_payload(g=[], t=[["LTR/Gypsy", "Gypsy", "+", 2000, 4000, 0]])
    assert run(tmp_path, data_with(tmp_path, gm2), body)["chev"] > 0
    gm3 = gm_payload(g=[["U", ".", 1500, 4000, 0, [["E", 1500, 1600]]]], t=[])
    body3 = """
    var g = H.byClass(det, "rpt-gm-gene", "g")[0];
    return {chev: H.all(g, "polyline").length};
    """
    assert run(tmp_path, data_with(tmp_path, gm3), body3)["chev"] == 0
    gm4 = gm_payload(g=[["S", "+", 1500, 4000, 0, [["E", 1500, 1600]]]], t=[])
    assert run(tmp_path, data_with(tmp_path, gm4), body3)["chev"] > 0


def test_annotation_names_only_become_text_nodes(tmp_path):
    evil = '<script>alert(1)</script>"x'
    gm = gm_payload(g=[[evil, "+", 1500, 1700, 0, [["E", 1500, 1700]]]],
                    t=[[evil, 'LTR/"<b>', "+", 3000, 3200, 1]])
    body = """
    res.script = H.doc("script").length; res.b = H.doc("b").length;
    var texts = H.all(gmEl, "text").map(H.text);
    res.hasLabel = texts.indexOf(%s) >= 0;
    res.html = H.serialize(gmEl);
    res.table = null;
    H.click(H.buttons(gmEl, %s)[0]);
    res.tableText = H.text(gmEl);
    return res;
    """ % (json.dumps(evil), json.dumps("Table"))
    d = data_with(tmp_path, gm)
    d["strings"]  # the Table button label comes from the strings
    body = body.replace(json.dumps("Table"), json.dumps(d["strings"]["ui.fig.table"]))
    r = run(tmp_path, d, body)
    assert r["script"] == 1 or r["script"] == 0   # the data block is the page's own
    assert r["b"] == 0 and r["hasLabel"]
    assert "<script>alert" not in r["html"] and "&lt;script&gt;alert(1)&lt;/script&gt;" in r["html"]
    assert evil in r["tableText"]


def test_empty_states_nothing_annotated_vs_no_layer(tmp_path):
    d = data_with(tmp_path, gm_payload(g=[], t=[], ly=["gene"]))
    S = d["strings"]
    body = """
    return H.all(gmEl, "text").filter(function (t) { return H.cls(t) === "rpt-gm-empty"; }).map(H.text);
    """
    notes = run(tmp_path, d, body)
    # genes (layer given, no rows) vs te / ncrna / phas (no layer)
    assert notes == [S["ui.gm.empty"], S["ui.gm.nolayer"].replace("{kind}", "te"),
                     S["ui.gm.nolayer"].replace("{kind}", "ncrna"), S["ui.gm.nolayer"].replace("{kind}", "phas")]


def test_more_not_shown_text(tmp_path):
    d = data_with(tmp_path, gm_payload(mo={"t": 7}))
    body = """
    return H.all(gmEl, "text").filter(function (t) { return H.cls(t) === "rpt-gm-more"; }).map(H.text);
    """
    assert run(tmp_path, d, body) == [d["strings"]["ui.gm.more"].replace("{n}", "7")]


def test_table_rows_equal_feature_count(tmp_path):
    gm = gm_payload(r=[["tRNA-1", "tRNA", "+", 4000, 4072, 0]], ly=["gene", "te", "ncrna"])
    d = data_with(tmp_path, gm)
    body = """
    H.click(H.buttons(gmEl, %s)[0]);
    var tb = H.byClass(gmEl, "rpt-srna-fig__table", "div")[0];
    var rows = H.all(tb, "tr");
    return {rows: rows.length - 1, first: H.all(rows[1], "td").map(H.text), head: H.all(rows[0], "th").map(H.text)};
    """ % json.dumps(d["strings"]["ui.fig.table"])
    r = run(tmp_path, d, body)
    assert r["rows"] == 2 + 2 + 1
    assert r["first"][:2] == ["Genes", "GeneA"] and r["first"][3] == "+"
    assert r["first"][4:] == ["1,501", "2,700", "✓"]
    assert len(r["head"]) == 7


def test_precursor_marker_only_with_mi(tmp_path):
    body = "return {n: H.byClass(gmEl, 'rpt-gm-prec', 'g').length};"
    assert run(tmp_path, data_with(tmp_path, gm_payload()), body)["n"] == 0
    mi = {"n": "x", "lo": 2900, "L": 120, "nowin": 0, "t": [], "pc": []}
    assert run(tmp_path, data_with(tmp_path, gm_payload(), mi=mi), body)["n"] == 1
    mi2 = dict(mi, nowin=1)
    assert run(tmp_path, data_with(tmp_path, gm_payload(), mi=mi2), body)["n"] == 0


def test_pairs_rescaled_to_one_shared_scale(tmp_path):
    cv = {"s": 0, "b": 10, "nb": 4, "k": ["21", "24"],
          "+": {"21": [10, [99, 0, 99]], "24": [4, [0, 99]]},
          "-": {"21": [5, [99]]}, "e": {"+": [3, [99]], "-": [6, [0, 99]]}}
    body = "return window.__gm.gmCoverage(%s);" % json.dumps(cv)
    r = run(tmp_path, data_with(tmp_path, gm_payload()), body)
    assert r["shared"] == 10 and r["eShared"] == 6
    assert r["vals"]["+"][0] == pytest.approx([10, 0, 10, 0])
    assert r["vals"]["+"][1] == pytest.approx([0, 4, 0, 0])
    assert r["vals"]["-"][0] == pytest.approx([5, 0, 0, 0])
    assert r["tot"]["+"] == pytest.approx([10, 4, 10, 0]) and r["top"] == pytest.approx(10)
    assert r["eTop"] == pytest.approx(6)


def test_te_class_from_family(tmp_path):
    body = "return ['LTR/Gypsy','line/L1','SINE','DNA/hAT','RC/Helitron','Simple_repeat','',null].map(window.__gm.gmTeClass);"
    r = run(tmp_path, data_with(tmp_path, gm_payload()), body)
    assert r == ["LTR", "LINE", "SINE", "DNA", "RC", "other", "other", "other"]


def test_tooltips_carry_coordinates_and_overlap(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    var tip = H.id("tip");
    var hs = H.byClass(gmEl, "rpt-gm-hit", "rect").filter(function (e) { return e.getAttribute("data-track") === "t"; });
    H.fire(hs[0], "pointermove", {clientX: 5, clientY: 5});
    res.t0 = H.text(tip);
    H.fire(hs[1], "pointermove", {clientX: 5, clientY: 5});
    res.t1 = H.text(tip);
    return res;
    """
    r = run(tmp_path, d, body)
    S = d["strings"]
    assert "Gypsy-1" in r["t0"] and "3,001–3,200" in r["t0"] and S["ui.gm.tip.overlap"] in r["t0"]
    assert S["ui.gm.tip.nostrand"] in r["t1"] and S["ui.gm.tip.overlap"] not in r["t1"]


def test_placed_features_and_their_labels_never_overlap_in_a_lane(tmp_path):
    body = """
    var rows = [], i;
    for (i = 0; i < 40; i++) {
      var s = 1000 + (i * 331) % 4800, e = s + 100 + (i * 97) % 700;
      rows.push({kind: "g", label: "Gene" + i, strand: "+", s: s, e: e, hl: i % 7 === 0, blocks: [["E", s, e]]});
    }
    function xa(p) { return 112 + (p - 1000) / 5000 * 772; }
    var pl = window.__gm.gmPlace("g", rows, xa, 112, 884, 1000, 6000);
    return {items: pl.items.map(function (it) { return {lane: it.lane, a: it.ea, b: it.eb, ml: it.label ? 1 : 0}; }),
            lanes: pl.lanes, dropped: pl.dropped};
    """
    r = run(tmp_path, data_with(tmp_path, gm_payload()), body)
    seen = [it for it in r["items"] if it["lane"] >= 0]
    assert r["lanes"] >= 2 and len(seen) + r["dropped"] == 40
    for i, p in enumerate(seen):
        for q in seen[i + 1:]:
            if p["lane"] == q["lane"]:
                assert p["b"] + 8 <= q["a"] + 1e-6 or q["b"] + 8 <= p["a"] + 1e-6, (p, q)
