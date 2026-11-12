"""Run report_assets/app.js under JavaScriptCore (osascript) with a stub DOM."""
import copy
import json
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_html_report_html import SUMMARY, data_of, full_sig, mk, zero_sig  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("osascript") is None,
                                reason="needs osascript (JavaScriptCore)")

HERE = os.path.dirname(os.path.abspath(__file__))
HARNESS = open(os.path.join(HERE, "jsc_harness.js"), encoding="utf-8").read()
APP = open(os.path.join(HERE, "..", "..", "report_assets", "app.js"), encoding="utf-8").read()
SVG_FIGS = "ABCDEFH"


class JsError(Exception):
    pass


def run_js(tmp_path, data, probe, app=APP, width=720):
    """Load app.js on the stub page, run `probe` (JS body using H, document), return its value."""
    cfg = {"app": app, "data": data, "probe": probe, "width": width}
    script = tmp_path / "run.js"
    script.write_text("var CFG = " + json.dumps(cfg) + ";\n" + HARNESS, encoding="utf-8")
    p = subprocess.run(["osascript", "-l", "JavaScript", str(script)],
                       capture_output=True, text=True, timeout=120)
    if p.returncode != 0:
        raise JsError(p.stderr.strip())
    out = json.loads(p.stdout.strip().splitlines()[-1])
    if not out["ok"]:
        raise JsError(out["error"])
    return out["result"]


# Overview state, then open the first clickable row; reports what rendered.
PROBE = """
var res = {};
res.tabs = H.all(H.id("tabs"), "button").map(H.text);
res.active0 = H.active();
res.ov = H.counts(H.id("sec-overview"));
res.ovText = H.text(H.id("sec-overview"));
res.navmeta = H.text(H.id("navmeta"));
var lst = H.id("sec-clusters");
res.list = H.counts(lst);
res.rows = H.all(lst, "tr").length - 1;
res.clickable = H.byClass(lst, "row", "tr").length;
res.listText = H.text(lst);
res.opened = false;
var rowEl = H.byClass(lst, "row", "tr")[0];
if (rowEl) {
  H.click(rowEl);
  var det = H.id("sec-detail");
  res.opened = true;
  res.active1 = H.active();
  res.detail = H.counts(det);
  // the checklist's row icons are not figures
  var ckEl = H.byClass(det, "rpt-srna-ck", "section")[0];
  res.detail.svg -= ckEl ? H.all(ckEl, "svg").length : 0;
  res.detailText = H.text(det);
  res.detailH = H.text(H.all(H.byClass(det, "rpt-srna-head")[0], "h2")[0]);
  res.strip = H.all(H.byClass(det, "rpt-srna-strip")[0], "div").filter(function (e) {
    return e.parent && H.cls(e.parent) === "rpt-srna-strip"; }).length;
  res.figs = H.byClass(det, "rpt-srna-fig", "section").map(function (f) { return f.getAttribute("data-fig"); });
  var seqs = H.byClass(det, "rpt-srna-seqs", "table");
  res.seqRows = seqs.length ? H.all(H.all(seqs[0], "tbody")[0], "tr").length : 0;
  res.notes = H.byClass(det, "nostrand", "text").map(H.text);
  H.click(H.buttons(det, STR_BACK)[0]);
  res.active2 = H.active();
}
res.script = H.doc("script").length;
res.bold = H.doc("b").length;
return res;
"""


def probe_for(data, body=PROBE):
    back = data["strings"]["ui.detail.back"]
    return "var STR_BACK = " + json.dumps(back) + ";\n" + body


def sigs(sel, kind="ok", make=None):
    out = {}
    for c in sel.detail:
        if make:
            out[c.name] = make()
        elif kind == "ok":
            out[c.name] = full_sig()
        elif kind == "none":
            out[c.name] = zero_sig()
        else:
            out[c.name] = {"failed": True}
    return out


def page(tmp_path, clusters, kind="ok", make=None, **kw):
    """Build the report with real-looking signals; return (sel, data json text, data dict)."""
    from report_data import select
    import report_html
    sel = select(clusters, kw.get("floor", 300.0), max_detail=kw.get("max_detail", 3000))
    out = tmp_path / "r.html"
    report_html.build_report(sel, kw.get("summary", SUMMARY), sigs(sel, kind, make), "t", str(out))
    html = out.read_text()
    d = data_of(html)
    return sel, json.dumps(d), d


def go(tmp_path, clusters, kind="ok", **kw):
    sel, text, d = page(tmp_path, clusters, kind, **kw)
    return sel, d, run_js(tmp_path, text, probe_for(d))


def check_common(d, r):
    S = d["strings"]
    assert r["tabs"] == [S["ui.about.title"], S["ui.tab.overview"], S["ui.tab.clusters"]]
    assert not any("Glossary" in t or "References" in t for t in r["tabs"])
    assert r["active0"] == ["overview"]
    assert r["navmeta"] == "t"
    assert r["script"] == 0 and r["bold"] == 0
    assert "glossary" not in d and "references" not in d


def check_page(r, figs=True):
    assert r["opened"] and r["active1"] == ["detail"] and r["active2"] == ["clusters"]
    assert r["strip"] >= 8                            # checks + library-context cells
    if figs:
        assert r["figs"] == list("ABCDEFGH")
        assert r["detail"]["svg"] == 7                # every figure but the sequences table
        assert 1 <= r["seqRows"] <= 10


def test_full_dataset(tmp_path):
    cl = [mk(f"c{i}", reads=400 + i, cls=("piRNA" if i % 3 else "siRNA"), review=i == 299)
          for i in range(300)]
    sel, d, r = go(tmp_path, cl, max_detail=3000)
    check_common(d, r)
    assert len(d["rows"]) == 300 and sum(x["d"] for x in d["rows"]) == 300
    assert r["rows"] == r["clickable"] == 50          # one page
    assert d["strings"]["ui.page.info"].format(p=1, pages=6, n=300) in r["listText"]
    assert r["ov"]["svg"] == 2                        # class bars + size histogram
    assert r["ov"]["path"] == 0                       # class bars use rects, no slices
    check_page(r)
    assert r["seqRows"] == 10
    assert r["detailH"] == "c299"                     # highest confidence tie -> most reads
    S = d["strings"]
    assert S["ui.review"] in r["detailText"]          # rv flag
    assert S["verdict.piRNA"] in r["detailText"]
    for gone in ("How to read", "Why not", "not needed for this class (", "Cautions"):
        assert gone not in r["detailText"]
    assert r["notes"] == []                           # both strands have reads


def test_paging_sorting_filtering(tmp_path):
    cl = [mk(f"c{i}", reads=400 + i) for i in range(120)]
    _, text, d = page(tmp_path, cl)
    probe = """
    var lst = H.id("sec-clusters");
    var n0 = H.all(lst, "tr").length - 1;
    H.click(H.buttons(lst, NEXT)[0]);
    var n1 = H.all(H.id("sec-clusters"), "tr").length - 1;
    H.click(H.buttons(H.id("sec-clusters"), NEXT)[0]);
    var n2 = H.all(H.id("sec-clusters"), "tr").length - 1;
    var sel = H.all(H.id("sec-clusters"), "select")[0];
    H.setValue(sel, "name");
    var n3 = H.all(H.id("sec-clusters"), "tr").length - 1;
    var inp = H.all(H.id("sec-clusters"), "input")[0];
    H.setValue(inp, 500);
    var n4 = H.all(H.id("sec-clusters"), "tr").length - 1;
    H.click(H.buttons(H.id("sec-clusters"), "piRNA (120)")[0]);
    return [n0, n1, n2, n3, n4];
    """
    nxt = json.dumps(d["strings"]["ui.page.next"])
    r = run_js(tmp_path, text, "var NEXT = " + nxt + ";\n" + probe)
    assert r[:3] == [50, 50, 20]
    assert r[3] == 50                                 # page reset after sort
    assert r[4] == 20                                 # floor 500 keeps reads 500..519


def test_empty_selection(tmp_path):
    sel, d, r = go(tmp_path, [mk(reads=10), mk("b", reads=20)], floor=300.0, summary={})
    check_common(d, r)
    assert d["rows"] == []
    assert r["ov"]["svg"] == 0                        # empty-state messages, no charts
    assert r["rows"] == 1 and r["clickable"] == 0     # the single "no rows" message row
    assert not r["opened"]


def test_capped_detail(tmp_path):
    cl = [mk(f"c{i}", reads=400 + i) for i in range(8)]
    sel, d, r = go(tmp_path, cl, max_detail=3)
    check_common(d, r)
    assert r["rows"] == 8 and r["clickable"] == 3
    assert d["meta"]["capped"] is True
    check_page(r)


def test_ambiguous_only(tmp_path):
    cl = [mk(f"a{i}", cls="ambiguous", gates={"piRNA": False, "miRNA": False},
             evidence="nearest=piRNA; size=0.6") for i in range(4)]
    sel, d, r = go(tmp_path, cl)
    check_common(d, r)
    assert r["rows"] == r["clickable"] == 4
    assert r["ov"]["svg"] == 2 and r["ov"]["circle"] == 0   # one class: one bar, no ring
    check_page(r)
    assert d["strings"]["ui.conf.amb"].format(cls="piRNA") in r["detailText"]
    assert d["strings"]["ui.ctx.conf.na"] in r["detailText"]


def test_plant_like(tmp_path):
    cl = [mk("h1", cls="hcsiRNA", gates={"hcsiRNA": True}),
          mk("h2", cls="hcsiRNA", gates={"hcsiRNA": True}),
          mk("p1", cls="phasiRNA", gates={"phasiRNA": True}, dicer_call=21),
          mk("m1", cls="miRNA", gates={"miRNA": True})]
    sel, d, r = go(tmp_path, cl)
    check_common(d, r)
    assert d["meta"]["kingdom"] == "plant"
    assert r["rows"] == 4 and r["ov"]["path"] == 0
    check_page(r)
    assert d["strings"]["ui.fig.pp.plant"] in r["detailText"]


def test_hostile_name(tmp_path):
    evil = '</script><b>"x"&'
    sel, text, d = page(tmp_path, [mk(evil, reads=900), mk("ok", reads=400)])
    r = run_js(tmp_path, text, probe_for(d, PROBE.replace("return res;", """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var det = H.id("sec-detail");
    var fig = H.byClass(det, "rpt-srna-fig", "section")[0];
    H.click(H.buttons(fig, DL)[0]);
    res.download = H.stats().downloads[0].name;
    res.blob = H.stats().lastBlob;
    res.script = H.doc("script").length;
    res.bold = H.doc("b").length;
    return res;""")).replace("var STR_BACK", "var DL = " + json.dumps(d["strings"]["ui.fig.svg"]) +
                          ";\nvar STR_BACK"))
    check_common(d, r)
    assert evil in r["listText"]
    assert r["detailH"] == evil
    assert r["script"] == 0 and r["bold"] == 0       # the name only ever became text
    assert r["download"] == "script_b_x_A.svg"
    assert "<b>" not in r["blob"] and "<script" not in r["blob"]
    check_page(r)


def test_old_tsv(tmp_path):
    cl = [mk(f"c{i}", e={}, pct_10a=None, has_evidence=False, confidence=None,
             size_purity=None, frac_top=None, pct_1u=None, pingpong_z=None,
             phase_pvalue=None) for i in range(3)]
    sel, d, r = go(tmp_path, cl)
    check_common(d, r)
    assert r["rows"] == 3
    assert r["ov"]["svg"] == 2 and "3" in r["ovText"]
    check_page(r)
    S = d["strings"]
    assert S["check.1u.title"] in r["detailText"] and S["status.unavailable"] in r["detailText"]


def test_strip_value_with_unavailable_level(tmp_path):
    # raw TSV value present, evidence score missing: value shown with a neutral marker
    cl = [mk("c0", e={}, has_evidence=False, pct_1u=0.81)]
    sel, text, d = page(tmp_path, cl)
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var st = H.byClass(H.id("sec-detail"), "rpt-srna-strip")[0];
    return st.children.map(function (c) {
      var lv = H.byClass(c, "lv")[0];
      return [H.text(c.children[0]), H.text(c.children[1]), lv ? H.cls(lv) : null];
    });""")
    one = [x for x in r if x[0] == d["strings"]["check.1u.title"]][0]
    assert one[1] == "81%" and one[2] == "lv unavailable"
    strand = [x for x in r if x[0] == d["strings"]["check.strand.title"]][0]
    assert strand[1] == d["strings"]["ui.unit.plus"].format(v=80)


def test_no_reads_cluster(tmp_path):
    sel, d, r = go(tmp_path, [mk("n1"), mk("n2")], kind="none")
    check_common(d, r)
    check_page(r, figs=False)
    assert r["detail"]["svg"] == 0 and r["figs"] == []
    assert d["strings"]["ui.fig.noreads"] in r["detailText"]
    assert "2" in r["ovText"]


def test_failed_region_cluster(tmp_path):
    sel, d, r = go(tmp_path, [mk("f1")], kind="failed")
    check_common(d, r)
    check_page(r, figs=False)
    assert r["detail"]["svg"] == 0
    assert d["strings"]["ui.fig.err"] in r["detailText"]
    assert d["strings"]["check.size.title"] in r["detailText"]
    assert r["ov"]["svg"] == 2


def test_tabs_switch(tmp_path):
    _, text, d = page(tmp_path, [mk()])
    r = run_js(tmp_path, text, """
    var b = H.all(H.id("tabs"), "button");
    var out = [];
    b.forEach(function (x) { H.click(x); out.push([H.active()[0], x.getAttribute("aria-selected")]); });
    out.push(H.byClass(H.id("tabs"), "rpt-srna-ind").length);
    return out;""")
    assert r[:3] == [["about", "true"], ["overview", "true"], ["clusters", "true"]]
    assert r[3] == 1


def test_about_tab(tmp_path):
    _, text, d = page(tmp_path, [mk()])
    r = run_js(tmp_path, text, """
    var a = H.id("sec-about");
    var links = H.all(a, "a");
    return {first: H.text(H.all(H.id("tabs"), "button")[0]), text: H.text(a),
      links: links.map(function (l) { return [l.getAttribute("href"), l.getAttribute("target"),
        l.getAttribute("rel")]; }), pre: H.text(H.all(a, "pre")[0]), svg: H.all(a, "svg").length};""")
    S, m = d["strings"], d["meta"]
    assert r["first"] == S["ui.about.title"]
    assert r["links"] == [[m["github"], "_blank", "noopener noreferrer"]]
    for t in (S["ui.about.name"], m["tagline"], S["ui.about.version"].format(v=m["version"]),
              S["ui.about.gh.title"], S["ui.about.gh.sub"]):
        assert t in r["text"]
    assert r["pre"] == "\n".join(m["install"])
    for gone in ("PyPI", "Paper", "Cite", "Citation", "DOI"):
        assert gone not in r["text"]


def test_row_keyboard_enter_and_space(tmp_path):
    _, text, d = page(tmp_path, [mk("a", reads=900), mk("b", reads=500)])
    r = run_js(tmp_path, text, """
    var row = H.byClass(H.id("sec-clusters"), "row", "tr")[0];
    var out = [row.getAttribute("tabindex")];
    H.key(row, "Tab");
    out.push(H.active()[0]);
    var ev = H.key(row, "Enter");
    out.push(H.active()[0]);
    H.click(H.buttons(H.id("sec-detail"), STR_BACK)[0]);
    row = H.byClass(H.id("sec-clusters"), "row", "tr")[1];
    ev = H.key(row, " ");
    out.push(H.active()[0], !!ev.prevented, H.text(H.all(H.id("sec-detail"), "h2")[0]));
    return out;""".replace("STR_BACK", json.dumps(d["strings"]["ui.detail.back"])))
    assert r == ["0", "overview", "detail", "detail", True, "b"]


# ---- figures: controls, tooltips, keyboard, strands, geometry ----

FIG_PROBE = """
H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
var det = H.id("sec-detail");
var out = {};
H.byClass(det, "rpt-srna-fig", "section").forEach(function (f) {
  var k = f.getAttribute("data-fig");
  var leg = H.byClass(f, "rpt-srna-legend", "ul");
  var hitEls = H.byClass(f, "rpt-srna-hit", "rect");
  var svg = H.all(f, "svg")[0];
  out[k] = {
    buttons: H.all(f, "button").map(H.text),
    legend: leg.length ? H.all(leg[0], "li").map(H.text) : [],
    title: H.text(H.byClass(f, "rpt-srna-fig__title")[0]),
    hits: hitEls.length,
    handlers: hitEls.every(function (e) {
      return ["pointermove", "pointerleave", "focus", "blur", "keydown"].every(function (t) { return H.has(e, t); });
    }),
    tab0: hitEls.filter(function (e) { return e.getAttribute("tabindex") === "0"; }).length,
    svg: svg ? H.serialize(svg) : null,
    width: svg ? svg.getAttribute("width") : null
  };
});
return out;
"""


def figs_of(tmp_path, clusters=None, make=None, width=720, **kw):
    clusters = clusters or [mk("c1", reads=900)]
    _, text, d = page(tmp_path, clusters, make=make, **kw)
    return d, run_js(tmp_path, text, FIG_PROBE, width=width)


def test_every_figure_has_controls_legend_and_targets(tmp_path):
    d, f = figs_of(tmp_path)
    S = d["strings"]
    assert sorted(f) == list("ABCDEFGH")
    nb = d["details"]["0"]["sig"]["nb"]
    want = {"A": 18, "B": 100, "C": nb, "D": 30, "E": nb, "F": 26}
    for k in SVG_FIGS:
        assert f[k]["buttons"] == [S["ui.fig.table"], S["ui.fig.svg"]], k
        assert f[k]["handlers"] and f[k]["hits"] > 0, k
        assert f[k]["width"] == "720", k
        if k in want:
            assert f[k]["hits"] == want[k], k
            assert f[k]["tab0"] == 1, k                # one roving tab stop per set
        assert re.search(r"^[A-H]", f[k]["title"]) and "n =" in f[k]["title"], k
    assert f["H"]["tab0"] == 2                         # two panels, two sets
    for k in "AF":
        assert f[k]["legend"][:4] == ["A", "C", "G", "U"], k
    assert S["ui.fig.period"].format(p=24) in f["B"]["legend"]   # piRNA: the dominant size
    assert f["G"]["buttons"] == [] and f["G"]["svg"] is None


def test_table_toggle_and_svg_download(tmp_path):
    _, text, d = page(tmp_path, [mk("c1", reads=900)])
    S = d["strings"]
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var out = {};
    H.byClass(H.id("sec-detail"), "rpt-srna-fig", "section").forEach(function (f) {
      var k = f.getAttribute("data-fig");
      var bs = H.all(f, "button");
      if (!bs.length) { return; }
      var plot = H.byClass(f, "rpt-srna-fig__plot")[0], tbl = H.byClass(f, "rpt-srna-fig__table")[0];
      H.click(bs[0]);
      var on = [bs[0].getAttribute("aria-pressed"), H.text(bs[0]), plot.hasAttribute("hidden"),
        tbl.hasAttribute("hidden"), H.text(H.all(tbl, "p")[0]), H.all(tbl, "tr").length,
        H.all(H.all(tbl, "tr")[1], "td").map(H.text)];
      H.click(bs[0]);
      var off = [bs[0].getAttribute("aria-pressed"), plot.hasAttribute("hidden"), tbl.hasAttribute("hidden")];
      H.click(bs[1]);
      var dl = H.stats().downloads[H.stats().downloads.length - 1];
      out[k] = {on: on, off: off, name: dl.name, blob: H.stats().lastBlob.slice(0, 60),
        body: H.all(document.body).length};
    });
    return out;""")
    assert sorted(r) == list(SVG_FIGS)
    for k, v in r.items():
        assert v["on"][:4] == ["true", S["ui.fig.chart"], True, False], k
        assert v["on"][4] in (S["ui.fig.precision"], S["ui.fig.exact"])
        assert v["on"][5] > 1
        assert v["off"] == ["false", False, True], k
        assert v["name"] == f"c1_{k}.svg"
        assert v["blob"].startswith('<svg xmlns="http://www.w3.org/2000/svg"')
        assert v["body"] == 1                          # the link was removed again
    # table values of A reproduce the payload within 1% of the panel maximum
    assert r["A"]["on"][6][0] == "18" and r["A"]["on"][6][1] == "+"


def test_size_table_matches_signal(tmp_path):
    sig = full_sig()
    _, text, d = page(tmp_path, [mk("c1", reads=900)], make=lambda: copy.deepcopy(sig))
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var f = H.byClass(H.id("sec-detail"), "rpt-srna-fig", "section")[0];
    H.click(H.all(f, "button")[0]);
    return H.all(H.byClass(f, "rpt-srna-fig__table")[0], "tr").slice(1).map(function (t) {
      return H.all(t, "td").map(H.text); });""")
    top = max(max(sig["sizes"][b][s]) for b in "ACGU" for s in "+-")
    assert len(r) == 36
    for row in r:
        i, s = int(row[0]) - 18, row[1]
        for j, b in enumerate("ACGU"):
            got = float(row[2 + j].replace(",", ""))
            assert abs(got - sig["sizes"][b][s][i]) <= 0.01 * top


def test_tooltips_on_pointer_focus_and_keys(tmp_path):
    _, text, d = page(tmp_path, [mk("c1", reads=900)])
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var tip = H.id("tip"), out = {};
    H.byClass(H.id("sec-detail"), "rpt-srna-fig", "section").forEach(function (f) {
      var k = f.getAttribute("data-fig");
      var hs = H.byClass(f, "rpt-srna-hit", "rect");
      if (!hs.length) { return; }
      var o = {};
      H.fire(hs[1], "pointermove", {clientX: 100, clientY: 100});
      o.move = [!tip.hasAttribute("hidden"), H.text(tip), H.cls(hs[1]), tip.style.transform];
      H.fire(hs[1], "pointerleave");
      o.leave = [tip.hasAttribute("hidden"), H.cls(hs[1])];
      hs[0].focus();
      o.focus = [!tip.hasAttribute("hidden"), hs[0].getAttribute("aria-label")];
      H.key(hs[0], "ArrowRight");
      o.right = [H.activeEl() === hs[1], hs[1].getAttribute("tabindex"), hs[0].getAttribute("tabindex")];
      H.key(hs[1], "End");
      o.end = H.activeEl() !== hs[1] && H.activeEl().getAttribute("tabindex") === "0";
      H.key(hs[hs.length - 1], "Escape");
      o.esc = tip.hasAttribute("hidden");
      out[k] = o;
    });
    return out;""")
    for k, o in r.items():
        assert o["move"][0] and o["move"][1] and "on" in o["move"][2].split(), k
        assert o["move"][3].startswith("translate("), k
        assert "NaN" not in o["move"][1] and "undefined" not in o["move"][1], k
        assert o["leave"] == [True, "rpt-srna-hit"], k
        assert o["focus"][0] and o["focus"][1], k
        assert o["right"] == [True, "0", "-1"], k
        assert o["end"] and o["esc"], k
    S = d["strings"]
    assert S["ui.strand.sense"] in r["B"]["move"][1] and S["ui.tip.dist"].format(x=2) in r["B"]["move"][1]


def one_strand(keep):
    def make():
        sig = full_sig()
        drop = "-" if keep == "+" else "+"
        for b in "ACGU":
            sig["sizes"][b][drop] = [0.0] * 18
        sig["phasing"][drop] = [0.0] * 100
        sig["cov"][drop] = [0.0] * len(sig["cov"][drop])
        sig["ends"][drop] = [0.0] * len(sig["ends"][drop])
        sig["logo"][drop] = [dict.fromkeys("ACGU", 0)] * 26
        sig["pingpong"] = [0.0] * 30
        sig["seqs"] = [q for q in sig["seqs"] if q[2] == keep][:3]
        return sig
    return make


@pytest.mark.parametrize("keep", ["+", "-"])
def test_single_strand_draws_an_empty_half(tmp_path, keep):
    d, f = figs_of(tmp_path, make=one_strand(keep))
    note = d["strings"]["ui.fig.nostrand"]
    for k in "ABCEF":
        svg = f[k]["svg"]
        assert svg.count(">" + note + "<") == 1, k
        assert "NaN" not in svg and "undefined" not in svg and "Infinity" not in svg, k
    assert d["strings"]["ui.fig.nopairs"] in f["D"]["svg"]
    # the note sits in the empty half
    root = ET.fromstring(f["A"]["svg"])
    zero = [e for e in root.iter() if e.get("class") == "zero"][0]
    y0 = float(zero.get("y1"))
    t = [e for e in root.iter() if e.text == note][0]
    assert (float(t.get("y")) > y0) == (keep == "+")


# ---- geometry, measured on the serialised SVG ----

NUM = re.compile(r"-?\d+(?:\.\d+)?")


def bbox(e):
    tag = e.tag.split("}")[-1]
    if tag == "rect":
        x, y, w, h = (float(e.get(k)) for k in ("x", "y", "width", "height"))
        return x, y, x + w, y + h
    if tag == "path":
        v = [float(x) for x in NUM.findall(e.get("d"))]
        xs, ys = v[0::2], v[1::2]
        return min(xs), min(ys), max(xs), max(ys)
    if tag == "line":
        x1, x2, y1, y2 = (float(e.get(k)) for k in ("x1", "x2", "y1", "y2"))
        return min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2)
    return None


def text_box(e):
    w = len(e.text or "") * 8.0
    x, y = float(e.get("x")), float(e.get("y"))
    a = e.get("text-anchor")
    x0 = x - w / 2 if a == "middle" else x - w if a == "end" else x
    return x0, y - 10, x0 + w, y + 3


def check_geometry(svg, W, arcs=False):
    root = ET.fromstring(svg)
    H = float(root.get("height"))
    assert float(root.get("width")) == W and root.get("viewBox") == f"0 0 {W} {H:g}"
    assert int(root.get("font-size")) >= 13
    for e in root.iter():
        assert e.get("font-size") is None or float(e.get("font-size")) >= 13
        b = None if arcs and e.tag.endswith("path") else bbox(e)   # arc flags are not points
        if b is not None:
            assert -0.01 <= b[0] and b[2] <= W + 0.01 and -0.01 <= b[1] and b[3] <= H + 0.01, (e.tag, b)
        if e.tag.endswith("text") and e.get("transform") is None:
            x0, y0, x1, y1 = text_box(e)
            assert x0 >= -1 and x1 <= W + 1 and y0 >= -1 and y1 <= H + 1, (e.text, x0, x1)
    # tick labels along each axis row never overlap
    rows = {}
    for e in root.iter():
        if "tk-x" in (e.get("class") or "").split():
            rows.setdefault(e.get("y"), []).append(text_box(e))
    for labs in rows.values():
        labs.sort()
        for a, b in zip(labs, labs[1:]):
            assert a[2] <= b[0], (a, b)
    marks = [e for e in root.iter() if any(c.startswith("mk-") for c in (e.get("class") or "").split())]
    return root, marks


@pytest.mark.parametrize("width", [720, 390, 330])
def test_geometry_inside_viewport_and_text_size(tmp_path, width):
    d, f = figs_of(tmp_path, width=width)
    for k in SVG_FIGS:
        root, marks = check_geometry(f[k]["svg"], width)
        assert marks, k
        for m in marks:
            if m.tag.endswith("path") and "mk-seg" in m.get("class"):
                x0, _, x1, _ = bbox(m)
                assert x1 - x0 <= 24.01                 # bars at most 24 px thick


def test_mirrored_halves_are_symmetric(tmp_path):
    def make():
        sig = full_sig()
        for b in "ACGU":
            sig["sizes"][b]["-"] = list(sig["sizes"][b]["+"])
        sig["cov"]["-"] = list(sig["cov"]["+"])
        sig["ends"]["-"] = list(sig["ends"]["+"])
        return sig
    d, f = figs_of(tmp_path, make=make)
    for k in "ACE":
        root, marks = check_geometry(f[k]["svg"], 720)
        zero = float([e for e in root.iter() if e.get("class") == "zero"][0].get("y1"))
        up = sorted((b[0], zero - b[3], zero - b[1])
                    for b in (bbox(m) for m in marks if "up" in m.get("class").split()))
        dn = sorted((b[0], b[1] - zero, b[3] - zero)
                    for b in (bbox(m) for m in marks if "dn" in m.get("class").split()))
        assert up and len(up) == len(dn), k
        for u, v in zip(up, dn):
            assert all(abs(a - b) <= 0.06 for a, b in zip(u, v)), (k, u, v)
        assert all(lo >= -0.01 for _, lo, _ in up), k  # up marks above the zero line


def test_redraw_on_resize(tmp_path):
    _, text, d = page(tmp_path, [mk("c1", reads=900)])
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var det = H.id("sec-detail");
    function figs() {   // not the 14 px checklist icons
      return H.all(det, "svg").map(function (s) { return s.getAttribute("width"); })
        .filter(function (w) { return w !== "14"; });
    }
    var before = figs();
    var n = H.observed();
    H.resize(410);
    var after = figs();
    return [before, n, after];""")
    assert set(r[0]) == {"720"} and r[1] == 9 and set(r[2]) == {"410"}   # 7 figures + 2 overview


# ---- fix round 1: numbers that could silently be wrong, lift, export, tabs ----

def fig_svg(tmp_path, clusters, make, letter, width=720):
    d, f = figs_of(tmp_path, clusters=clusters, make=make, width=width)
    return d, ET.fromstring(f[letter]["svg"])


def plot_x(i, n_minus_1, W=720):
    return LEFT + i / n_minus_1 * (W - RIGHT - LEFT)


LEFT, RIGHT = 64, 14


def test_phasing_offset_and_period_positions(tmp_path):
    def make():
        sig = full_sig()
        for s in "+-":
            sig["phasing"][s] = [0.0] * 100
            sig["phasing"][s][20] = 100.0             # a peak at offset 21
        return sig
    cl = [mk("p1", cls="phasiRNA", gates={"phasiRNA": True}, dicer_call=21)]
    d, root = fig_svg(tmp_path, cl, make, "B")
    line = [e for e in root.iter() if "mk-line" in (e.get("class") or "")][0]
    pts = [tuple(map(float, p.split(","))) for p in re.findall(r"[ML]([\d.]+,[\d.]+)", line.get("d"))]
    assert len(pts) == 100
    peak = min(pts, key=lambda p: p[1])
    assert abs(peak[0] - plot_x(20, 99)) < 0.02        # offset 21 is index 20
    xs = sorted({round(float(e.get("x1")), 2) for e in root.iter() if e.get("class") == "period"})
    assert xs == [round(plot_x(k - 1, 99), 2) for k in (21, 42, 63, 84)]


def test_pingpong_red_bar_is_ten_and_plant_is_grey(tmp_path):
    def make():
        sig = full_sig()
        sig["pingpong"] = [50.0] * 30
        return sig
    _, root = fig_svg(tmp_path, [mk("c1", reads=900)], make, "D")
    bars = {e.get("data-i"): e.get("fill") for e in root.iter() if "mk-bar" in (e.get("class") or "")}
    assert len(bars) == 30
    assert [k for k, v in bars.items() if v == "#e41a1c"] == ["9"]
    assert {v for k, v in bars.items() if k != "9"} == {"#999999"}
    plant = [mk("h1", cls="hcsiRNA", gates={"hcsiRNA": True})]
    _, root = fig_svg(tmp_path, plant, make, "D")
    fills = {e.get("fill") for e in root.iter() if "mk-bar" in (e.get("class") or "")}
    assert fills == {"#999999"}


def test_logo_puts_the_largest_base_on_top(tmp_path):
    def make():
        sig = full_sig()
        for s in "+-":
            sig["logo"][s] = list(sig["logo"][s])
            sig["logo"][s][0] = {"A": 0.1, "C": 0.6, "G": 0.25, "U": 0.05}
        return sig
    _, root = fig_svg(tmp_path, [mk("c1", reads=900)], make, "F")
    segs = [e for e in root.iter() if "mk-seg" in (e.get("class") or "") and e.get("data-i") == "0"]
    top = sorted(segs, key=lambda e: bbox(e)[1])[:4]  # the sense panel's column, top down
    assert [e.get("fill") for e in top] == ["#3B54A2", "#F9A31B", "#108041", "#EA2627"]


def test_coverage_coordinates_are_one_based(tmp_path):
    _, text, d = page(tmp_path, [mk("c1", reads=900)])
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var out = {};
    H.byClass(H.id("sec-detail"), "rpt-srna-fig", "section").forEach(function (f) {
      var k = f.getAttribute("data-fig");
      if (k !== "C" && k !== "E") { return; }
      H.fire(H.byClass(f, "rpt-srna-hit", "rect")[0], "pointermove", {clientX: 5, clientY: 5});
      var head = H.text(H.byClass(H.id("tip"), "hd")[0]);
      H.click(H.all(f, "button")[0]);
      out[k] = [head, H.text(H.all(H.all(H.byClass(f, "rpt-srna-fig__table")[0], "tbody")[0], "td")[0])];
    });
    return out;""")
    s0, b = d["details"]["0"]["sig"]["cv"]["s"], d["details"]["0"]["sig"]["cv"]["b"]
    want = d["strings"]["ui.tip.pos"].format(start=f"{s0 + 1:,}", end=f"{s0 + b:,}")
    assert r["C"] == [want, want] and r["E"] == [want, want]


def test_coverage_ticks_are_one_based_bases(tmp_path):
    d, root = fig_svg(tmp_path, [mk("c1", reads=900)], None, "C")
    cv = d["details"]["0"]["sig"]["cv"]
    s0, end = cv["s"], cv["s"] + d["details"]["0"]["sig"]["nb"] * cv["b"]

    def xa(p0):                                         # x of a 0-based position
        return LEFT + (p0 - s0) / (end - s0) * (720 - RIGHT - LEFT)
    ticks = [(int(e.text.replace(",", "")), float(e.get("x"))) for e in root.iter()
             if "tk-x" in (e.get("class") or "") and e.get("text-anchor") == "middle"]
    assert len(ticks) >= 2
    step = ticks[1][0] - ticks[0][0]
    for p, x in ticks:
        assert p > s0 and p % step == 0                 # round 1-based positions
        assert abs(x - xa(p - 1)) < 0.02, (p, x)        # at the left edge of base p


def test_logo_title_shows_exact_nd_only(tmp_path):
    import math
    d, f = figs_of(tmp_path)
    nd = d["details"]["0"]["sig"]["nd"]
    want = d["strings"]["ui.fig.n"].format(n=f"{int(math.floor(nd + 0.5)):,}")
    assert want in f["F"]["title"] and f["F"]["title"].count("n\u00a0=") == 1
    assert "n\u00a0=" not in f["F"]["svg"]                # no per-strand n in the panels


def test_hovered_column_lifts_and_clears(tmp_path):
    _, text, d = page(tmp_path, [mk("c1", reads=900)])
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var figs = {};
    H.byClass(H.id("sec-detail"), "rpt-srna-fig", "section").forEach(function (f) {
      figs[f.getAttribute("data-fig")] = f; });
    function state(f) {
      var svg = H.all(f, "svg")[0];
      var marks = H.all(svg).filter(function (e) { return e.getAttribute && e.getAttribute("data-i") !== null; });
      return [svg.classList.contains("lift"),
        marks.filter(function (e) { return e.classList.contains("cur"); }).map(function (e) {
          return e.getAttribute("data-i"); }).filter(function (v, i, a) { return a.indexOf(v) === i; })];
    }
    var out = {};
    "ADFH".split("").forEach(function (k) {
      var hs = H.byClass(figs[k], "rpt-srna-hit", "rect"), o = {};
      var i = k === "D" ? 9 : 2;
      H.fire(hs[i], "pointermove", {clientX: 9, clientY: 9});
      o.on = state(figs[k]);
      H.fire(hs[i], "pointerleave", {pointerType: "mouse"});
      o.off = state(figs[k]);
      hs[i].focus();
      o.focus = state(figs[k]);
      H.key(hs[i], "Escape");
      o.esc = state(figs[k]);
      H.fire(hs[i], "pointermove", {clientX: 9, clientY: 9});
      H.fire(hs[i], "pointerleave", {pointerType: "touch"});
      o.touch = !H.id("tip").hasAttribute("hidden");
      out[k] = o;
    });
    var c = H.byClass(figs.C, "rpt-srna-hit", "rect");
    H.fire(c[3], "pointermove", {clientX: 9, clientY: 9});
    out.C = [state(figs.C)[0], H.cls(c[3])];
    return out;""")
    for k, i in (("A", "2"), ("D", "9"), ("F", "2"), ("H", "0:2")):
        o = r[k]
        if k == "H" and not o["on"][0]:
            continue                                   # an empty histogram bin has no bar
        assert o["on"] == [True, [i]], k
        assert o["off"] == [False, []] and o["esc"] == [False, []], k
        assert o["focus"] == [True, [i]], k
        assert o["touch"], k                           # a lifted finger keeps the tooltip
    assert r["C"] == [False, "rpt-srna-hit on"]          # area figures: the wash only


def test_svg_export_has_no_hit_layer(tmp_path):
    _, text, d = page(tmp_path, [mk("c1", reads=900)])
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var f = H.byClass(H.id("sec-detail"), "rpt-srna-fig", "section")[0];
    var hs = H.byClass(f, "rpt-srna-hit", "rect");
    H.fire(hs[2], "pointermove", {clientX: 9, clientY: 9});
    H.click(H.all(f, "button")[1]);
    return [H.stats().lastBlob, hs.length, H.byClass(f, "rpt-srna-hit", "rect").length,
      H.all(f, "svg")[0].classList.contains("lift")];""")
    blob = r[0]
    for bad in ("rpt-srna-hit", "tabindex", "data-h", 'role="img"', "lift", " cur", '"cur'):
        assert bad not in blob, bad
    assert "mk-seg" in blob and r[1] == r[2] == 18 and r[3]   # the page keeps its layer


def test_tablist_keyboard(tmp_path):
    _, text, d = page(tmp_path, [mk()])
    r = run_js(tmp_path, text, """
    var b = H.all(H.id("tabs"), "button"), out = [];
    function st() {
      return [H.active()[0], H.activeEl() ? H.activeEl().getAttribute("id") : null,
        b.map(function (x) { return x.getAttribute("tabindex"); }).join("")];
    }
    out.push(b.map(function (x) { return [x.getAttribute("id"), x.getAttribute("aria-controls")]; }));
    out.push(["about", "overview", "clusters", "detail"].map(function (s) {
      var e = H.id("sec-" + s); return [e.getAttribute("role"), e.getAttribute("aria-labelledby")]; }));
    out.push(st());
    var ev = H.key(b[1], "ArrowRight"); out.push(st().concat([!!ev.prevented]));
    H.key(b[2], "ArrowRight"); out.push(st());
    H.key(b[0], "ArrowLeft"); out.push(st());
    H.key(b[2], "Home"); out.push(st());
    H.key(b[0], "End"); out.push(st());
    return out;""")
    assert r[0] == [["tab-about", "sec-about"], ["tab-overview", "sec-overview"],
                    ["tab-clusters", "sec-clusters"]]
    assert r[1] == [["tabpanel", "tab-about"], ["tabpanel", "tab-overview"],
                    ["tabpanel", "tab-clusters"], ["tabpanel", "tab-clusters"]]
    assert r[2][0] == "overview" and r[2][2] == "-10-1"
    assert r[3] == ["clusters", "tab-clusters", "-1-10", True]
    assert r[4][:2] == ["about", "tab-about"]          # wraps
    assert r[5][:2] == ["clusters", "tab-clusters"]
    assert r[6][:2] == ["about", "tab-about"] and r[7][:2] == ["clusters", "tab-clusters"]


def test_ambiguous_strip_is_neutral_but_shows_values(tmp_path):
    cl = [mk("a0", cls="ambiguous", gates={"piRNA": False, "miRNA": False},
             evidence="nearest=piRNA; size=0.6")]
    _, text, d = page(tmp_path, cl)
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var st = H.byClass(H.id("sec-detail"), "rpt-srna-strip")[0];
    return st.children.map(function (c) {
      var lv = H.byClass(c, "lv")[0];
      return [H.text(c.children[1]), lv ? H.cls(lv) : null];
    });""")
    S = d["strings"]
    checks = [x for x in r if x[1] is not None]
    assert len(checks) == 7
    assert {x[1] for x in checks} == {"lv context"}     # nothing defines "ambiguous"
    vals = [x[0] for x in checks]
    assert vals[:4] == [S["ui.unit.nt"].format(v=26), "80%", "50%", S["ui.unit.plus"].format(v=80)]
    assert S["ui.na"] not in vals


def test_focus_survives_a_redraw(tmp_path):
    _, text, d = page(tmp_path, [mk("c1", reads=900)])
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    var f = H.byClass(H.id("sec-detail"), "rpt-srna-fig", "section")[0];
    H.byClass(f, "rpt-srna-hit", "rect")[3].focus();
    H.resize(500);
    var a = H.activeEl(), svg = H.all(f, "svg")[0];
    return [a.getAttribute("data-h"), a.getAttribute("tabindex"), a.parent === svg, svg.getAttribute("width")];""")
    assert r == ["0:3", "0", True, "500"]


def test_overview_charts_in_real_pixels(tmp_path):
    _, text, d = page(tmp_path, [mk(f"c{i}", reads=400 + i, cls=("piRNA" if i % 2 else "siRNA"))
                                 for i in range(6)])
    for w in (720, 330):
        r = run_js(tmp_path, text, """
        return H.all(H.id("sec-overview"), "svg").map(function (s) { return H.serialize(s); });""", width=w)
        assert len(r) == 2
        for svg in r:
            check_geometry(svg, w, arcs=True)


def test_inner_html_is_refused_by_the_stub(tmp_path):
    _, text, _ = page(tmp_path, [mk()])
    with pytest.raises(JsError, match="innerHTML"):
        run_js(tmp_path, text, 'document.createElement("div").innerHTML = "<b>x</b>"; return 1;')


# ---- the harness itself must fail loudly ----
def test_harness_detects_thrown_error(tmp_path):
    _, text, _ = page(tmp_path, [mk()])
    with pytest.raises(JsError, match="boom"):
        run_js(tmp_path, text, "return 1;", app=APP.replace('"use strict";',
               '"use strict"; throw new Error("boom");', 1))


def test_harness_detects_syntax_error(tmp_path):
    _, text, _ = page(tmp_path, [mk()])
    with pytest.raises(JsError, match="SyntaxError"):
        run_js(tmp_path, text, "return 1;", app=APP + "\n)(")


def test_harness_detects_missing_dom_method(tmp_path):
    _, text, _ = page(tmp_path, [mk()])
    bad = APP.replace('"use strict";', '"use strict"; document.querySelectorAll("x");', 1)
    with pytest.raises(JsError, match="TypeError"):
        run_js(tmp_path, text, "return 1;", app=bad)


def test_harness_detects_probe_failure(tmp_path):
    _, text, _ = page(tmp_path, [mk()])
    with pytest.raises(JsError):
        run_js(tmp_path, text, 'H.click(H.id("sec-about")); return 1;')


CK_PROBE = """
H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
var det = H.id("sec-detail");
var ck = H.byClass(det, "rpt-srna-ck", "section")[0];
if (!ck) { return null; }
var order = det.children.map(function (c) { return H.cls(c); });
return {text: H.text(ck), head: H.text(H.all(ck, "h3")[0]),
        rows: (H.all(ck, "tbody")[0] ? H.all(H.all(ck, "tbody")[0], "tr") : []).map(function (tr) {
          return H.all(tr, "td").map(H.text).concat(H.text(H.all(tr, "th")[0])); }),
        svg: H.all(ck, "svg").length,
        before: H.byClass(det, "rpt-srna-figs").length ? true : false};
"""


def ck_page(tmp_path, clusters):
    _, text, d = page(tmp_path, clusters)
    return d, run_js(tmp_path, text, CK_PROBE)


def test_checklist_rows_icons_and_text(tmp_path):
    d, r = ck_page(tmp_path, [mk("a", cls="piRNA", dicer_call=27, start=0, end=5000)])
    S = d["strings"]
    assert r["head"] == S["ck.title"].format(cls="piRNA")
    assert len(r["rows"]) == 4
    assert r["rows"][0][-1] == S["check.size.title"]
    assert S["ck.st.1"] in r["rows"][0][2]            # result is text, not only colour
    assert r["svg"] == 4                              # and an icon per row
    assert S["ck.thr.length"].format(n=1000) in r["text"]
    assert S["ck.note"] in r["text"]


def test_checklist_failed_rows_say_not_met(tmp_path):
    d, r = ck_page(tmp_path, [mk("a", cls="ambiguous", dicer_call=19, evidence="nearest=piRNA; x=1",
                                  e={k: 0.0 for k in ("1u", "5a", "fractop_bid", "pingpong", "phase")})])
    S = d["strings"]
    assert r["head"] == S["ck.title.amb"].format(cls="piRNA")
    assert [x[2] for x in r["rows"]] == [S["ck.st.0"], S["ck.st.0"], S["ck.st.3"], S["ck.st.3"]]


def test_checklist_mirna_has_one_line(tmp_path):
    d, r = ck_page(tmp_path, [mk("m", cls="miRNA")])
    assert d["strings"]["ck.mirna"] in r["text"] and r["rows"] == []


def test_checklist_plant_and_old_tsv(tmp_path):
    d, r = ck_page(tmp_path, [mk("h", cls="hcsiRNA", dicer_call=24, gates={"hcsiRNA": True},
                                  start=0, end=300)])
    assert len(r["rows"]) == 2
    cl = [mk("o", cls="siRNA", dicer_call=21, e={}, has_evidence=False)]
    d, r = ck_page(tmp_path, cl)
    assert r["rows"][1][2] == d["strings"]["ck.st.2"]


def test_checklist_sits_between_strip_and_figures(tmp_path):
    _, text, d = page(tmp_path, [mk("a", cls="piRNA", dicer_call=27, start=0, end=5000)])
    r = run_js(tmp_path, text, """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    return H.id("sec-detail").children.map(function (c) { return H.cls(c); });""")
    i = r.index("rpt-srna-ck")
    assert r[i - 1] == "rpt-srna-strip" and r[i + 1] == "rpt-srna-figs"
