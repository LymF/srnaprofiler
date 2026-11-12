"""Run report_assets/app.js under JavaScriptCore (osascript) with a stub DOM."""
import json
import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
from test_html_report_html import SUMMARY, build, data_of, full_sig, mk, zero_sig  # noqa: E402

pytestmark = pytest.mark.skipif(shutil.which("osascript") is None,
                                reason="needs osascript (JavaScriptCore)")

HERE = os.path.dirname(os.path.abspath(__file__))
HARNESS = open(os.path.join(HERE, "jsc_harness.js"), encoding="utf-8").read()
APP = open(os.path.join(HERE, "..", "..", "report_assets", "app.js"), encoding="utf-8").read()


class JsError(Exception):
    pass


def run_js(tmp_path, data, probe, app=APP):
    """Load app.js on the stub page, run `probe` (JS body using H, document), return its value."""
    cfg = {"app": app, "data": data, "probe": probe}
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
res.gloss = H.byClass(H.id("sec-glossary"), "rpt-card").length;
res.refs = H.all(H.id("sec-references"), "li").length;
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
  res.detailText = H.text(det);
  res.detailH = H.byClass(det, "rpt-section-title")[0].textContent;
  H.click(H.buttons(det, document_back())[0]);
  res.active2 = H.active();
}
function document_back() { return STR_BACK; }
return res;
"""


def probe_for(data):
    back = data["strings"]["ui.detail.back"]
    return "var STR_BACK = " + json.dumps(back) + ";\n" + PROBE


def sigs(sel, kind="ok"):
    out = {}
    for c in sel.detail:
        if kind == "ok":
            out[c.name] = full_sig()
        elif kind == "none":
            out[c.name] = zero_sig()
        else:
            out[c.name] = {"failed": True}
    return out


XFAIL = pytest.mark.xfail(strict=True, reason="front end lands in Task 5")


def page(tmp_path, clusters, kind="ok", **kw):
    """Build the report with real-looking signals; return (sel, data json text, data dict)."""
    from report_data import select
    import report_html
    sel = select(clusters, kw.get("floor", 300.0), max_detail=kw.get("max_detail", 3000))
    out = tmp_path / "r.html"
    report_html.build_report(sel, kw.get("summary", SUMMARY), sigs(sel, kind), "t", str(out))
    html = out.read_text()
    d = data_of(html)
    return sel, json.dumps(d), d


def go(tmp_path, clusters, kind="ok", **kw):
    sel, text, d = page(tmp_path, clusters, kind, **kw)
    return sel, d, run_js(tmp_path, text, probe_for(d))


def check_common(d, r):
    assert len(r["tabs"]) == 4
    assert r["active0"] == ["overview"]
    assert r["refs"] == len(d["references"])
    assert r["gloss"] == len(d["glossary"])
    assert r["navmeta"] == "t"


@XFAIL
def test_full_dataset(tmp_path):
    cl = [mk(f"c{i}", reads=400 + i, cls=("piRNA" if i % 3 else "siRNA")) for i in range(300)]
    sel, d, r = go(tmp_path, cl, max_detail=3000)
    check_common(d, r)
    assert len(d["rows"]) == 300 and sum(x["d"] for x in d["rows"]) == 300
    assert r["rows"] == r["clickable"] == 50          # one page
    assert d["strings"]["ui.page.info"].format(p=1, pages=6, n=300) in r["listText"]
    assert r["ov"]["svg"] == 2                        # donut + size histogram
    assert r["ov"]["path"] == 2                       # two classes in the donut
    assert r["opened"] and r["active1"] == ["detail"] and r["active2"] == ["clusters"]
    assert r["detail"]["svg"] == 4 and r["detail"]["rect"] > 0
    assert r["detailH"] == "c299"                     # highest confidence tie -> most reads


@XFAIL
def test_paging_sorting_filtering(tmp_path):
    cl = [mk(f"c{i}", reads=400 + i) for i in range(120)]
    _, text, d = page(tmp_path, cl)
    probe = """
    var lst = H.id("sec-clusters");
    var n0 = H.all(lst, "tr").length - 1;
    H.click(H.buttons(lst, NEXT)[0]);
    lst = H.id("sec-clusters");
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


@XFAIL
def test_empty_selection(tmp_path):
    sel, d, r = go(tmp_path, [mk(reads=10), mk("b", reads=20)], floor=300.0, summary={})
    check_common(d, r)
    assert d["rows"] == []
    assert r["ov"]["svg"] == 0                        # empty-state messages, no charts
    assert r["rows"] == 1 and r["clickable"] == 0     # the single "no rows" message row
    assert not r["opened"]


@XFAIL
def test_capped_detail(tmp_path):
    cl = [mk(f"c{i}", reads=400 + i) for i in range(8)]
    sel, d, r = go(tmp_path, cl, max_detail=3)
    check_common(d, r)
    assert r["rows"] == 8 and r["clickable"] == 3
    assert d["meta"]["capped"] is True
    assert r["opened"] and r["detail"]["svg"] == 4


@XFAIL
def test_ambiguous_only(tmp_path):
    cl = [mk(f"a{i}", cls="ambiguous", gates={"piRNA": False, "miRNA": False}) for i in range(4)]
    sel, d, r = go(tmp_path, cl)
    check_common(d, r)
    assert r["rows"] == r["clickable"] == 4
    assert r["ov"]["svg"] == 2 and r["ov"]["circle"] == 1   # one class: full ring
    assert r["opened"] and r["detail"]["svg"] == 4


@XFAIL
def test_plant_like(tmp_path):
    cl = [mk("h1", cls="hcsiRNA", gates={"hcsiRNA": True}),
          mk("h2", cls="hcsiRNA", gates={"hcsiRNA": True}),
          mk("p1", cls="phasiRNA", gates={"phasiRNA": True}, dicer_call=21),
          mk("m1", cls="miRNA", gates={"miRNA": True})]
    sel, d, r = go(tmp_path, cl)
    check_common(d, r)
    assert d["meta"]["kingdom"] == "plant"
    assert r["rows"] == 4 and r["ov"]["path"] == 3
    assert r["opened"] and r["detail"]["svg"] == 4


@XFAIL
def test_hostile_name(tmp_path):
    evil = '</script><b>"x"&'
    sel, d, r = go(tmp_path, [mk(evil, reads=900), mk("ok", reads=400)])
    check_common(d, r)
    assert evil in r["listText"]
    assert r["detailH"] == evil
    assert r["detail"]["svg"] == 4


@XFAIL
def test_old_tsv(tmp_path):
    cl = [mk(f"c{i}", e={}, pct_10a=None, has_evidence=False, confidence=None,
             size_purity=None, frac_top=None, pct_1u=None, pingpong_z=None,
             phase_pvalue=None) for i in range(3)]
    sel, d, r = go(tmp_path, cl)
    check_common(d, r)
    assert r["rows"] == 3
    assert r["opened"] and r["detail"]["svg"] == 4


@XFAIL
def test_no_reads_cluster(tmp_path):
    sel, d, r = go(tmp_path, [mk("n1"), mk("n2")], kind="none")
    check_common(d, r)
    assert r["opened"] and r["detail"]["svg"] == 0
    assert d["strings"]["ui.fig.noreads"] in r["detailText"]


@XFAIL
def test_failed_region_cluster(tmp_path):
    sel, d, r = go(tmp_path, [mk("f1")], kind="failed")
    check_common(d, r)
    assert r["opened"] and r["detail"]["svg"] == 0
    assert d["strings"]["ui.fig.err"] in r["detailText"]


@XFAIL
def test_tabs_switch(tmp_path):
    _, text, d = page(tmp_path, [mk()])
    r = run_js(tmp_path, text, """
    var b = H.all(H.id("tabs"), "button");
    var out = [];
    b.forEach(function (x) { H.click(x); out.push(H.active()[0]); });
    return out;""")
    assert r == ["overview", "clusters", "glossary", "references"]


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
        run_js(tmp_path, text, 'H.click(H.id("sec-glossary")); return 1;')
