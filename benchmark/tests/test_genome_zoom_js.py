"""Task 6: zoom, pan and locus focus in the genomic-context figure (stub DOM, JavaScriptCore)."""
import json
import os
import re
import shutil
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_html_report_js import APP, run_js  # noqa: E402
import test_genome_figure_js as t5  # noqa: E402
from test_genome_figure_js import data_with, gm_payload  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("osascript") is None,
                                reason="needs osascript (JavaScriptCore)")

HOOK = ("window.__gm = {gmZoom: gmZoom, gmPan: gmPan, gmLocusView: gmLocusView, gmBinsInView: gmBinsInView};"
        "\n  buildTabs();\n  renderAbout();")
APPX = APP.replace("  buildTabs();\n  renderAbout();", HOOK, 1)
APP_FILE = os.environ.get("GM_APP")
if APP_FILE:
    APPX = open(APP_FILE, encoding="utf-8").read()
assert APPX != APP or APP_FILE

WS, WE = t5.WS, t5.WE      # 1000, 6000
PRE = """
var plot = H.byClass(gmEl, "rpt-srna-fig__plot", "div")[0];
var rng = H.byClass(gmEl, "rpt-srna-fig__n", "span")[0];
var S = %s;
function btn(key) { return H.all(gmEl, "button").filter(function (b) { return b.getAttribute("aria-label") === S[key]; })[0]; }
function dis(key) { return btn(key).hasAttribute("disabled"); }
function view() {
  var m = /:([\\d,]+)-([\\d,]+),/.exec(H.text(rng));
  return [parseInt(m[1].replace(/,/g, ""), 10) - 1, parseInt(m[2].replace(/,/g, ""), 10)];
}
function go() { H.frames(); }
function svgOf() { return H.all(plot, "svg")[0]; }
function genes() { return H.byClass(svgOf(), "rpt-gm-gene", "g"); }
function geneY(g) { return +H.all(g, "line")[0].getAttribute("y1"); }
// x of the first and last plot pixel at the page width of the probe
var X0 = 112, X1 = 884;
function px(frac) { return 10 + X0 + frac * (X1 - X0); }   // clientX; the stub plot box starts at 10
"""


def run(tmp_path, d, body="return res;", width=900):
    pre = PRE % json.dumps(d["strings"])
    return t5.run_js(tmp_path, json.dumps(d), t5.OPEN + pre + body, app=APPX, width=width)


def helpers(tmp_path, expr):
    d = data_with(tmp_path, gm_payload())
    return run(tmp_path, d, "return " + expr + ";")


# ---------------------------------------------------------------- pure helpers
def test_zoom_in_halves_the_span_around_the_anchor(tmp_path):
    r = helpers(tmp_path, """{
      mid: window.__gm.gmZoom([1000, 6000], [1000, 6000], 2, 0.5),
      left: window.__gm.gmZoom([1000, 6000], [1000, 6000], 2, 0),
      right: window.__gm.gmZoom([1000, 6000], [1000, 6000], 2, 1),
      q: window.__gm.gmZoom([2000, 4000], [1000, 6000], 2, 0.25),
      out: window.__gm.gmZoom([2000, 4000], [1000, 6000], 0.5, 0.5)}""")
    assert r["mid"] == [2250, 4750] and r["left"] == [1000, 3500] and r["right"] == [3500, 6000]
    assert r["q"] == [2250, 3250]            # the position at 1/4 of the view stays at 1/4
    assert r["out"] == [1000, 5000]          # span 4000 around 3000, shifted inside the window


def test_zoom_is_clamped_to_40_bp_and_to_the_window(tmp_path):
    r = helpers(tmp_path, """{
      min: window.__gm.gmZoom([3000, 3050], [1000, 6000], 2, 0.5),
      min2: window.__gm.gmZoom([3000, 3040], [1000, 6000], 8, 0.1),
      big: window.__gm.gmZoom([1000, 6000], [1000, 6000], 0.5, 0.5),
      big2: window.__gm.gmZoom([2000, 4000], [1000, 6000], 0.001, 0.9),
      short: window.__gm.gmZoom([500, 520], [500, 520], 2, 0.5),
      short2: window.__gm.gmZoom([500, 520], [500, 520], 0.5, 0.5),
      flat: window.__gm.gmZoom([500, 500], [500, 500], 2, 0.5),
      noanchor: window.__gm.gmZoom([1000, 6000], [1000, 6000], 2)}""")
    assert r["min"] == [3005, 3045] and r["min2"][1] - r["min2"][0] == 40
    assert r["big"] == [1000, 6000] and r["big2"] == [1000, 6000]
    assert r["short"] == [500, 520] and r["short2"] == [500, 520]      # window < 40 bp: stays whole
    assert r["flat"] == [500, 500] and r["noanchor"] == [2250, 4750]


def test_zoom_keeps_the_pointer_position_fixed(tmp_path):
    body = """
    var out = [], v = [1000, 6000], W = [1000, 6000], fr = [0, 0.1, 0.37, 0.5, 0.93, 1], i, f;
    for (i = 0; i < fr.length; i++) {
      for (f = 0; f < 4; f++) {
        var pos = v[0] + fr[i] * (v[1] - v[0]);
        var n = window.__gm.gmZoom(v, W, 2, fr[i]);
        var pos2 = n[0] + fr[i] * (n[1] - n[0]);
        out.push([pos, pos2, n[0] >= 1000, n[1] <= 6000]);
        v = n;
      }
      v = [1000, 6000];
    }
    return out;
    """
    for pos, pos2, lo, hi in helpers_body(tmp_path, body):
        assert abs(pos - pos2) <= 1 and lo and hi


def helpers_body(tmp_path, body):
    d = data_with(tmp_path, gm_payload())
    return run(tmp_path, d, body)


def test_pan_moves_and_clamps_at_the_edges(tmp_path):
    r = helpers(tmp_path, """{
      a: window.__gm.gmPan([2000, 3000], [1000, 6000], 200),
      left: window.__gm.gmPan([1100, 2100], [1000, 6000], -500),
      right: window.__gm.gmPan([4900, 5900], [1000, 6000], 5000),
      whole: window.__gm.gmPan([1000, 6000], [1000, 6000], 300),
      frac: window.__gm.gmPan([2000, 3000], [1000, 6000], 12.4)}""")
    assert r["a"] == [2200, 3200] and r["left"] == [1000, 2000] and r["right"] == [5000, 6000]
    assert r["whole"] == [1000, 6000] and r["frac"] == [2012, 3012]


def test_locus_view_is_locus_plus_ten_percent_clamped(tmp_path):
    r = helpers(tmp_path, """{
      mid: window.__gm.gmLocusView([1000, 6000], [3000, 3100]),
      edge: window.__gm.gmLocusView([1000, 6000], [1000, 1100]),
      tail: window.__gm.gmLocusView([1000, 6000], [5950, 6000]),
      tiny: window.__gm.gmLocusView([1000, 6000], [3000, 3010]),
      huge: window.__gm.gmLocusView([1000, 6000], [500, 9000]),
      short: window.__gm.gmLocusView([500, 520], [505, 510])}""")
    assert r["mid"] == [2990, 3110] and r["edge"] == [1000, 1110] and r["tail"] == [5945, 6000]
    assert r["tiny"] == [2985, 3025]          # a 10 bp locus gets the 40 bp minimum, centred
    assert r["huge"] == [1000, 6000] and r["short"] == [500, 520]


def test_bins_in_view(tmp_path):
    cv = {"s": 1000, "b": 50, "nb": 100}
    body = "var cv = %s, B = window.__gm.gmBinsInView; return {full: B(cv, [1000, 6000]), one: B(cv, [1000, 1050])," \
           " part: B(cv, [1025, 1075]), edge: B(cv, [1050, 1100]), end: B(cv, [5990, 6000]), out: B(cv, [7000, 8000])," \
           " before: B(cv, [0, 900]), empty: B({s: 0, b: 10, nb: 0}, [0, 100])};" % json.dumps(cv)
    r = helpers_body(tmp_path, body)
    assert r["full"] == {"first": 0, "last": 99} and r["one"] == {"first": 0, "last": 0}
    assert r["part"] == {"first": 0, "last": 1} and r["edge"] == {"first": 1, "last": 1}
    assert r["end"] == {"first": 99, "last": 99}
    for k in ("out", "before", "empty"):
        assert r[k]["first"] > r[k]["last"], k


# ---------------------------------------------------------------- controls and state
def test_buttons_are_real_buttons_with_aria_labels(tmp_path):
    d = data_with(tmp_path, gm_payload())
    S = d["strings"]
    body = """
    var keys = ["ui.gm.zoom.in", "ui.gm.zoom.out", "ui.gm.zoom.reset", "ui.gm.pan.left", "ui.gm.pan.right", "ui.gm.locus"];
    var tools = H.byClass(gmEl, "rpt-srna-fig__tools", "div")[0];
    res.labels = keys.map(function (k) { var b = btn(k); return b ? [b.tagName, b.getAttribute("aria-label"), b.getAttribute("type")] : null; });
    res.inTools = keys.map(function (k) { return H.all(tools, "button").indexOf(btn(k)) >= 0; });
    res.order = H.all(tools, "button").length;
    res.range = H.text(rng);
    res.plotTab = plot.getAttribute("tabindex");
    return res;
    """
    r = run(tmp_path, d, body)
    keys = ["ui.gm.zoom.in", "ui.gm.zoom.out", "ui.gm.zoom.reset", "ui.gm.pan.left", "ui.gm.pan.right", "ui.gm.locus"]
    assert r["labels"] == [["button", S[k], "button"] for k in keys]
    assert all(r["inTools"]) and r["order"] == 8      # six controls, Table, Download SVG
    assert r["range"] == S["ui.gm.range"].replace("{chrom}", "c1").replace("{start}", "1,001") \
        .replace("{end}", "6,000").replace("{n}", "5,000") or "1,001-6,000, 5,000" in r["range"].replace(" ", " ")
    assert r["plotTab"] == "0"


def test_disabled_states_at_the_limits(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    var K = ["ui.gm.zoom.in", "ui.gm.zoom.out", "ui.gm.zoom.reset", "ui.gm.pan.left", "ui.gm.pan.right", "ui.gm.locus"];
    function st() { return K.map(dis); }
    res.full = st();
    var steps = 0;
    while (!dis("ui.gm.zoom.in") && steps < 30) { H.click(btn("ui.gm.zoom.in")); steps++; }
    go();
    res.steps = steps; res.min = st(); res.minView = view();
    H.click(btn("ui.gm.pan.left")); go();       // the centre zoom is not at an edge: both pans are free
    res.mid = st();
    var n = 0;
    while (!dis("ui.gm.pan.left") && n < 400) { H.click(btn("ui.gm.pan.left")); n++; }
    go();
    res.atLeft = st(); res.leftView = view();
    n = 0;
    while (!dis("ui.gm.pan.right") && n < 800) { H.click(btn("ui.gm.pan.right")); n++; }
    go();
    res.atRight = st(); res.rightView = view();
    H.click(btn("ui.gm.zoom.reset")); go();
    res.reset = st(); res.resetView = view();
    return res;
    """
    r = run(tmp_path, d, body)
    # order: in, out, reset, left, right, locus
    assert r["full"] == [False, True, True, True, True, False]
    assert r["min"][0] is True and r["min"][1:3] == [False, False] and r["minView"][1] - r["minView"][0] == 40
    assert r["mid"][3:5] == [False, False]
    assert r["atLeft"][3] is True and r["atLeft"][4] is False and r["leftView"][0] == WS
    assert r["atRight"][4] is True and r["atRight"][3] is False and r["rightView"][1] == WE
    assert r["reset"] == r["full"] and r["resetView"] == [WS, WE]


def test_locus_button_and_range_text(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    H.click(btn("ui.gm.locus"));
    res.pending = H.pendingFrames();
    res.afterClick = H.text(rng);     // the text does not wait for the frame
    go();
    res.view = view(); res.pendingAfter = H.pendingFrames();
    res.locusDis = dis("ui.gm.locus");
    res.svgRange = H.all(svgOf(), "text").map(H.text);
    return res;
    """
    r = run(tmp_path, d, body)
    assert r["view"] == [2990, 3110] and r["pending"] == 1 and r["pendingAfter"] == 0
    assert "2,991-3,110, 120" in r["afterClick"].replace(" ", " ") and r["locusDis"] is True
    assert "3,000" in r["svgRange"] or "3,050" in r["svgRange"]     # the ruler follows the view


def test_redraws_are_throttled_to_one_frame(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    var before = svgOf();
    H.click(btn("ui.gm.zoom.in")); H.click(btn("ui.gm.zoom.in")); H.click(btn("ui.gm.pan.right"));
    res.same = svgOf() === before;
    res.pending = H.pendingFrames();
    go();
    res.changed = svgOf() !== before;
    res.view = view();
    return res;
    """
    r = run(tmp_path, d, body)
    assert r["same"] and r["pending"] == 1 and r["changed"]
    assert r["view"][1] - r["view"][0] == 1250


# ---------------------------------------------------------------- keyboard, wheel, drag
def test_keyboard_when_the_plot_has_focus(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    plot.focus();
    res.k = [];
    function key(k) { var ev = H.key(plot, k); go(); res.k.push([k, view(), !!ev.prevented]); }
    key("+");  key("=");  key("ArrowRight");  key("ArrowLeft");  key("ArrowLeft");  key("-");  key("0");
    key("-");  key("ArrowLeft");  key("x");
    return res;
    """
    r = run(tmp_path, d, body)
    k = {i: v for i, v in enumerate(r["k"])}
    assert k[0][1] == [2250, 4750] and k[0][2]
    assert k[1][1] == [2875, 4125]
    span = k[1][1][1] - k[1][1][0]
    assert k[2][1][0] == k[1][1][0] + round(span * 0.2) and k[2][2]        # 20 % of the view
    assert k[3][1] == k[1][1] and k[4][1][0] == k[1][1][0] - round(span * 0.2)
    assert k[5][1][1] - k[5][1][0] == 2 * span                                # zoom out
    assert k[6][1] == [WS, WE] and k[6][2]                                    # 0 resets
    assert k[7][1] == [WS, WE] and k[8][1] == [WS, WE]                        # at the limits: nothing moves
    assert k[9][1] == [WS, WE] and not k[9][2]                                # other keys are left alone


def test_arrow_keys_on_a_feature_still_walk_the_features(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    var hit = H.byClass(svgOf(), "rpt-gm-hit", "rect")[0];
    H.fire(hit, "keydown", {key: "ArrowRight"}); go();
    res.view = view();
    return res;
    """
    assert run(tmp_path, d, body)["view"] == [WS, WE]


def test_wheel_zooms_around_the_pointer_only_when_focused_or_pinching(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    var ev0 = H.fire(plot, "wheel", {clientX: px(0.25), clientY: 40, deltaY: -100}); go();
    res.unfocused = [view(), !!ev0.prevented];
    var ev1 = H.fire(plot, "wheel", {clientX: px(0.25), clientY: 40, deltaY: -100, ctrlKey: true}); go();
    res.ctrl = [view(), !!ev1.prevented];
    H.click(btn("ui.gm.zoom.reset")); go();
    plot.focus();
    var v0 = view(), pos0 = v0[0] + 0.25 * (v0[1] - v0[0]);
    var ev2 = H.fire(plot, "wheel", {clientX: px(0.25), clientY: 40, deltaY: -100}); go();
    var v1 = view(); res.focused = [v1, !!ev2.prevented, pos0, v1[0] + 0.25 * (v1[1] - v1[0])];
    var ev3 = H.fire(plot, "wheel", {clientX: px(0.9), clientY: 40, deltaY: 100}); go();
    res.out = [view(), !!ev3.prevented];
    for (var i = 0; i < 8; i++) { H.fire(plot, "wheel", {clientX: px(1), clientY: 40, deltaY: 100}); }
    go(); res.outAll = view();
    return res;
    """
    r = run(tmp_path, d, body)
    assert r["unfocused"] == [[WS, WE], False]                    # the page keeps scrolling
    assert r["ctrl"][1] and r["ctrl"][0][1] - r["ctrl"][0][0] == 2500       # a pinch (ctrl + wheel) always zooms
    v, prevented, pos0, pos1 = r["focused"]
    assert prevented and v[1] - v[0] == 2500 and abs(pos0 - pos1) <= 1
    assert r["out"][1] and r["out"][0][1] - r["out"][0][0] == 5000
    assert r["outAll"] == [WS, WE]


def test_drag_pans_and_a_click_does_not(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    H.click(btn("ui.gm.zoom.in")); H.click(btn("ui.gm.zoom.in")); go();
    var v0 = view();
    H.fire(plot, "pointerdown", {clientX: px(0.5), clientY: 50, button: 0, pointerId: 1});
    H.fire(plot, "pointermove", {clientX: px(0.5) + 1, clientY: 50, pointerId: 1});
    go(); res.tiny = view();
    H.fire(plot, "pointermove", {clientX: px(0.5) + (X1 - X0) / 4, clientY: 50, pointerId: 1}); go();
    res.drag = [v0, view()];
    H.fire(plot, "pointerup", {clientX: px(0.75), clientY: 50, pointerId: 1});
    H.fire(plot, "pointermove", {clientX: px(0.1), clientY: 50, pointerId: 1}); go();
    res.after = view();
    H.fire(plot, "pointerdown", {clientX: px(0.5), clientY: 50, button: 0, pointerId: 1});
    H.fire(plot, "pointermove", {clientX: px(0.5) + 4000, clientY: 50, pointerId: 1}); go();
    res.clamped = view();
    H.fire(plot, "pointerup", {pointerId: 1});
    return res;
    """
    r = run(tmp_path, d, body)
    v0 = r["drag"][0]
    span = v0[1] - v0[0]
    assert r["tiny"] == v0                                      # below the drag threshold
    assert abs(r["drag"][1][0] - (v0[0] - span / 4)) <= 1          # dragging right moves the view left
    assert r["after"] == r["drag"][1]                            # pointer up ends the drag
    assert r["clamped"][0] == WS


def test_view_range_and_lanes_survive_zoom(tmp_path):
    gm = gm_payload(L=[3200, 3250],
                    g=[["GeneA", "+", 3000, 3120, 1, [["E", 3000, 3120]]],
                       ["Gb", "-", 3100, 3300, 0, [["E", 3100, 3300]]],
                       ["Gc", "+", 3050, 3400, 0, [["E", 3050, 3400]]]], t=[])
    d = data_with(tmp_path, gm)
    body = """
    var full = genes().map(geneY);
    H.click(btn("ui.gm.zoom.in")); H.click(btn("ui.gm.zoom.in")); go();
    res.n2 = genes().length;
    res.two = genes().map(geneY);
    H.click(btn("ui.gm.zoom.reset")); go();
    H.click(btn("ui.gm.locus")); go();
    res.locus = genes().map(geneY);
    res.full = full;
    res.h = [H.all(svgOf(), "svg").length, +svgOf().getAttribute("height")];
    return res;
    """
    r = run(tmp_path, d, body)
    assert len(set(r["full"])) == 3                  # three overlapping genes, three lanes
    # zoomed twice around the centre (3500): view [2875, 4125] still shows all three, in the same lanes
    assert r["n2"] == 3 and r["two"] == r["full"]
    assert r["locus"] == r["full"][1:]               # locus view [3195, 3255]: Gb and Gc, same lanes


def test_a_feature_outside_the_view_is_not_drawn_and_one_cut_by_the_edge_gets_a_marker(tmp_path):
    gm = gm_payload(g=[], t=[["LTR/Gypsy", "Gypsy-1", "+", 2000, 3050, 0],
                              ["DNA/hAT", "hAT-2", ".", 4500, 4600, 0]], L=[3000, 3040])
    d = data_with(tmp_path, gm)
    body = """
    H.click(btn("ui.gm.locus")); go();
    var bars = H.byClass(svgOf(), "rpt-gm-bar", "g");
    res.bars = bars.length;
    res.clip = bars.map(function (b) { return b.getAttribute("data-clip"); });
    res.cont = bars.length ? H.all(bars[0], "path").length : 0;
    res.view = view();
    return res;
    """
    r = run(tmp_path, d, body)
    assert r["bars"] == 1 and r["clip"] == ["lr"] and r["cont"] == 2     # cut at both edges, hAT is outside
    assert r["view"] == [2996, 3044]


def test_coverage_is_clipped_and_rescaled_to_the_bins_in_view(tmp_path):
    d = data_with(tmp_path, gm_payload())
    body = """
    function bars() {
      var out = [], tot = 0;
      H.byClass(svgOf(), "rpt-gm-cov", "path").forEach(function (p) {
        var dd = p.getAttribute("d"), m = dd.match(/M/g);
        tot += m ? m.length : 0;
        (dd.match(/M[^Z]*Z/g) || []).forEach(function (seg) {
          var ys = seg.match(/,(-?[\\d.]+)/g).map(function (s) { return +s.slice(1); });
          out.push(Math.max.apply(null, ys) - Math.min.apply(null, ys));
        });
      });
      return {count: tot, tallest: Math.max.apply(null, out.concat([0]))};
    }
    res.full = bars();
    // zoom (anchor at the left edge) to 40 bp, then step right until only the third bin is in view
    plot.focus();
    for (var i = 0; i < 9; i++) { H.fire(plot, "wheel", {clientX: px(0), clientY: 40, deltaY: -100}); }
    go(); res.min = view();
    for (i = 0; i < 13; i++) { H.key(plot, "ArrowRight"); }
    go(); res.bin2 = view(); res.zoomed = bars();
    res.edges = H.byClass(svgOf(), "rpt-gm-end", "path").length;
    return res;
    """
    r = run(tmp_path, d, body)
    assert r["min"] == [1000, 1040]
    assert r["bin2"][0] >= 1100 and r["bin2"][1] <= 1150           # inside the third 50 bp bin only
    assert r["zoomed"]["count"] == 3                                 # + 21 nt, + 24 nt, - 21 nt of that bin
    # the third bin is 1.8 of 14 at full scale; rescaled it is the tallest thing on the axis
    assert r["full"]["tallest"] > r["zoomed"]["tallest"] * 0 and r["zoomed"]["tallest"] > 25
    assert r["edges"] == 0                                           # no 5' ends in that bin


def test_table_lists_only_features_overlapping_the_view(tmp_path):
    gm = gm_payload(L=[3200, 3250], r=[["tRNA-1", "tRNA", "+", 4000, 4072, 0]], ly=["gene", "te", "ncrna"],
                    g=[["GeneA", "+", 3000, 3120, 1, [["E", 3000, 3120]]], ["Gb", "-", 3100, 3300, 0, [["E", 3100, 3300]]]],
                    t=[["LTR/Gypsy", "Gypsy-1", "+", 3000, 3200, 1], ["DNA/hAT", "hAT-2", ".", 3300, 3500, 0]])
    d = data_with(tmp_path, gm)
    S = d["strings"]
    body = """
    var tb = H.byClass(gmEl, "rpt-srna-fig__table", "div")[0];
    var tbtn = H.buttons(gmEl, %s)[0];
    function names() { return H.all(tb, "tr").slice(1).map(function (r) { return H.all(r, "td").map(H.text)[1]; }); }
    H.click(tbtn); res.full = names();
    H.click(tbtn);                       // back to the chart
    H.click(btn("ui.gm.locus"));         // the view changes while the table is closed (frame not yet run)
    H.click(tbtn); res.locus = names();
    H.click(btn("ui.gm.zoom.in")); go();   // and while it is open
    res.open = names(); res.view = view();
    return res;
    """ % json.dumps(S["ui.fig.table"])
    r = run(tmp_path, d, body)
    assert r["full"] == ["GeneA", "Gb", "LTR/Gypsy", "DNA/hAT", "tRNA-1"]
    # locus view [3195, 3255]: Gb (3100-3300) and Gypsy-1 (3000-3200) overlap it, GeneA (ends 3120) does not
    assert r["locus"] == ["Gb", "LTR/Gypsy"]
    assert r["view"][1] - r["view"][0] == 40 and r["open"] == ["Gb"]      # 3215-3255 is only inside Gb


def test_download_svg_is_the_current_view(tmp_path):
    gm = gm_payload(g=[["GeneA", "+", 3000, 3120, 1, [["E", 3000, 3120]]], ["Gb", "-", 4100, 4300, 0, [["E", 4100, 4300]]]], t=[])
    d = data_with(tmp_path, gm)
    body = """
    var dl = H.buttons(gmEl, %s)[0];
    function count() { H.click(dl); return (H.stats().lastBlob.match(/rpt-gm-gene/g) || []).length; }
    res.full = count();
    H.click(btn("ui.gm.locus"));         // no frame has run: the download must still be the new view
    res.locus = count();
    res.name = H.stats().downloads.length;
    return res;
    """ % json.dumps(d["strings"]["ui.fig.svg"])
    r = run(tmp_path, d, body)
    assert r["full"] == 2 and r["locus"] == 1 and r["name"] == 2     # locus view: only GeneA


def test_names_with_markup_stay_text_after_zoom(tmp_path):
    evil = '<script>alert(1)</script>"x'
    gm = gm_payload(L=[3000, 3100], g=[[evil, "+", 2990, 3110, 0, [["E", 2990, 3110]]]],
                    t=[[evil, 'LTR/"<b>', "+", 3000, 3200, 1]])
    d = data_with(tmp_path, gm)
    body = """
    H.click(btn("ui.gm.locus")); H.click(btn("ui.gm.zoom.in")); go();
    res.script = H.doc("script").length; res.b = H.doc("b").length;
    res.html = H.serialize(gmEl);
    res.has = H.all(svgOf(), "text").map(H.text).indexOf(%s) >= 0 || H.all(svgOf(), "rect").length > 0;
    return res;
    """ % json.dumps(evil)
    r = run(tmp_path, d, body)
    assert r["b"] == 0 and "<script>alert" not in r["html"] and "&lt;script&gt;" in r["html"]


def test_plain_report_without_gm_is_unchanged(tmp_path):
    d = data_with(tmp_path)
    r = t5.run_js(tmp_path, json.dumps(d), t5.OPEN + "return {has: res.has, n: res.figs.length};", app=APPX, width=900)
    assert r["has"] is False and r["n"] >= 1


def test_zoom_buttons_labels_are_not_in_the_unused_allowlist():
    src = open(os.path.join(os.path.dirname(__file__), "..", "..", "report_assets", "app.js"), encoding="utf-8").read()
    for k in ("zoom.in", "zoom.out", "zoom.reset", "pan.left", "pan.right"):
        assert re.search(r'"ui\.gm\.%s"' % re.escape(k), src), k
