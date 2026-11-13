"""Annotation evidence in the HTML report: data layer, payload keys, legend, page."""
import json
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import report_html  # noqa: E402
from report_data import load_clusters  # noqa: E402
from test_html_report_html import build, data_of, mk  # noqa: E402

TEXT = "te:LTR/Gypsy:antisense:0.82|gene:exon:sense:0.10"
LEGEND = "never changes the class"
HEAD = ["name", "chrom", "start", "end", "reads", "reads_rpm", "cls", "confidence",
        "dicer_call", "size_purity", "frac_top", "pct_1u", "pct_10a", "pingpong_z",
        "phase_pvalue", "multimap_frac", "unique_frac"]


def write_tsv(path, annotated):
    cols = HEAD + (["annotation", "annotation_support"] if annotated else []) + ["review"]
    vals = ["c1", "chr1", "1000", "1100", "500", "5", "piRNA", "0.9", "26", "0.8", "0.8",
            "0.8", "0.5", "3", "0.5", "0.1", "0.9"]
    if annotated:
        vals += [TEXT, "conflicts"]
    vals += ["True"]
    path.write_text("\t".join(cols) + "\n" + "\t".join(vals) + "\n")
    return str(path)


def test_load_clusters_reads_annotation_columns(tmp_path):
    c = load_clusters(write_tsv(tmp_path / "a.tsv", True))[0]
    assert c.annotation == TEXT and c.annotation_support == "conflicts" and c.review
    c = load_clusters(write_tsv(tmp_path / "p.tsv", False))[0]
    assert c.annotation == "" and c.annotation_support == "none"


def test_row_has_an_only_when_annotated(tmp_path):
    sel, html, _ = build(tmp_path, [mk("a", annotation=TEXT, annotation_support="conflicts"),
                                    mk("b"), mk("c", annotation=TEXT, annotation_support="none")])
    rows = {r["n"]: r for r in data_of(html)["rows"]}
    assert rows["a"]["an"] == "conflicts"
    assert "an" not in rows["b"] and "an" not in rows["c"]


def test_detail_carries_annotation_text(tmp_path):
    _, html, _ = build(tmp_path, [mk("a", annotation=TEXT, annotation_support="supports"),
                                  mk("b", reads=400)])
    d = data_of(html)
    idx = {r["n"]: str(i) for i, r in enumerate(d["rows"])}
    assert d["details"][idx["a"]]["an"] == TEXT
    assert "an" not in d["details"][idx["b"]]


def test_report_html_mentions_annotation_legend_only_when_used(tmp_path):
    _, html, _ = build(tmp_path, [mk("a", annotation=TEXT, annotation_support="supports")])
    d = data_of(html)
    assert d["meta"]["an"] and LEGEND in d["strings"]["ui.an.legend"]
    _, plain, _ = build(tmp_path, [mk("a")])
    d = data_of(plain)
    assert "an" not in d["meta"] and "an_sum" not in d["meta"]
    assert all("an" not in r for r in d["rows"])
    assert all("an" not in x for x in d["details"].values())


@pytest.mark.skipif(shutil.which("osascript") is None, reason="needs osascript")
def test_page_shows_text_chip_and_legend(tmp_path):
    from test_html_report_js import go, run_js, probe_for, page
    _, text, d = page(tmp_path, [mk("a", annotation=TEXT, annotation_support="conflicts")])
    probe = """
    var lst = H.id("sec-clusters"); var res = {};
    res.list = H.text(lst);
    H.click(H.byClass(lst, "row", "tr")[0]);
    res.detail = H.text(H.id("sec-detail"));
    res.chips = H.byClass(H.id("sec-detail"), "rpt-srna-an", "span").map(H.text);
    return res;"""
    r = run_js(tmp_path, text, probe_for(d, probe))
    assert TEXT in r["detail"] and r["chips"] == ["conflicts"] and LEGEND in r["detail"]
    assert r["chips"] == ["conflicts"]
    _, text, d = page(tmp_path, [mk("a")])
    r = run_js(tmp_path, text, probe_for(d, probe))
    assert LEGEND not in r["detail"] and r["chips"] == []


def test_an_sum_counts_shown_loci_per_class(tmp_path):
    cl = [mk("a", annotation_support="supports"), mk("b", annotation_support="conflicts"),
          mk("c", annotation_support="supports"),
          mk("d", cls="siRNA", annotation_support="neutral"),
          mk("hidden", reads=5, annotation_support="conflicts")]  # below the floor
    sel, html, _ = build(tmp_path, cl)
    assert len(sel.shown) == 4
    assert data_of(html)["meta"]["an_sum"] == {"piRNA": [2, 1, 0], "siRNA": [0, 0, 1]}


def test_row_annotation_support_is_whitelisted(tmp_path):
    cl = [mk("a", annotation_support='x" onclick="y'), mk("b", annotation_support="supports")]
    _, html, _ = build(tmp_path, cl)
    d = data_of(html)
    rows = {r["n"]: r for r in d["rows"]}
    assert "an" not in rows["a"] and rows["b"]["an"] == "supports"
    assert d["meta"]["an_sum"] == {"piRNA": [1, 0, 0]}


def test_an_sum_absent_without_annotation(tmp_path):
    _, html, _ = build(tmp_path, [mk("a"), mk("b", annotation_support="none")])
    assert "an_sum" not in data_of(html)["meta"]


@pytest.mark.skipif(shutil.which("osascript") is None, reason="needs osascript")
def test_overview_renders_summary_table_with_text_nodes(tmp_path):
    from test_html_report_js import run_js, probe_for, page
    cl = [mk("a", annotation_support="supports"), mk("b", annotation_support="conflicts"),
          mk("d", cls="siRNA", annotation_support="neutral")]
    _, text, d = page(tmp_path, cl)
    probe = """
    var ov = H.id("sec-overview");
    return {text: H.text(ov), tables: H.all(ov, "table").length,
            ths: H.all(ov, "th").map(H.text), tds: H.all(ov, "td").map(H.text)};"""
    r = run_js(tmp_path, text, probe_for(d, probe))
    S = d["strings"]
    assert S["ui.an.title"] in r["text"]
    for k in ("ui.an.th.class", "ui.an.th.supports", "ui.an.th.conflicts", "ui.an.th.neutral"):
        assert S[k] in r["ths"]
    assert r["tds"] == ["piRNA", "1", "1", "0", "siRNA", "0", "0", "1"]
    _, text, d = page(tmp_path, [mk("a")])
    r = run_js(tmp_path, text, probe_for(d, probe))
    assert S["ui.an.title"] not in r["text"] and r["tds"] == []


@pytest.mark.skipif(shutil.which("osascript") is None, reason="needs osascript")
def test_detail_chip_class_never_takes_unknown_value(tmp_path):
    from test_html_report_js import run_js, probe_for, page
    _, text, d = page(tmp_path, [mk("a", annotation=TEXT, annotation_support="supports")])
    d["rows"][0]["an"] = 'x" onclick="y'
    import json as _j
    probe = """
    H.click(H.byClass(H.id("sec-clusters"), "row", "tr")[0]);
    return H.byClass(H.id("sec-detail"), "rpt-srna-an", "span").map(function (e) { return e.className; });"""
    r = run_js(tmp_path, _j.dumps(d), probe_for(d, probe))
    assert r == []


LIST_PROBE = """
var lst = H.id("sec-clusters");
function snap() {
  var rows = H.byClass(lst, "row", "tr");
  return {names: rows.map(function (t) { return H.text(H.all(H.all(t, "td")[0], "strong")[0]); }),
          chips: H.byClass(lst, "rpt-srna-an", "span").map(function (e) { return e.className; }),
          ths: H.all(lst, "th").map(H.text),
          cells: rows.map(function (t) { return H.all(t, "td").length; }),
          shown: H.text(lst)};
}
function ctl(k) { return H.all(lst, "select").filter(function (e) { return e.getAttribute("data-ctl") === k; })[0] || null; }
var res = {first: snap(), hasAn: !!ctl("an")};
if (ctl("an")) {
  res.opts = H.all(ctl("an"), "option").map(function (o) { return [o.getAttribute("value"), H.text(o)]; });
  H.setValue(ctl("an"), "conflicts");
  res.conf = snap();
  res.focus = H.activeEl() && H.activeEl().getAttribute("data-ctl");
  H.click(H.all(lst, "button").filter(function (b) { return H.text(b).indexOf("piRNA") === 0; })[0]);
  res.both = snap();
  H.click(H.all(lst, "button").filter(function (b) { return H.text(b).indexOf("siRNA") === 0; })[0]);
  res.sirna = snap();
}
return res;"""

AN_ROWS = [("a", "piRNA", "supports"), ("b", "piRNA", "conflicts"), ("c", "siRNA", "conflicts"),
           ("d", "siRNA", "neutral"), ("e", "piRNA", "none")]


@pytest.mark.skipif(shutil.which("osascript") is None, reason="needs osascript")
def test_list_annotation_column_and_filter(tmp_path):
    from test_html_report_js import run_js, probe_for, page
    cl = [mk(n, cls=c, annotation_support=s, reads=500 - i)
          for i, (n, c, s) in enumerate(AN_ROWS)]
    _, text, d = page(tmp_path, cl)
    S = d["strings"]
    r = run_js(tmp_path, text, probe_for(d, LIST_PROBE))
    assert r["hasAn"] and S["ui.an.filter"] in r["first"]["shown"]
    assert S["ui.an.filter"] in r["first"]["ths"]
    assert [o[0] for o in r["opts"]] == ["all", "supports", "conflicts", "neutral"]
    assert sorted(r["first"]["chips"]) == ["rpt-srna-an conflicts", "rpt-srna-an conflicts",
                                           "rpt-srna-an neutral", "rpt-srna-an supports"]
    assert set(r["first"]["cells"]) == {7}            # row without a verdict keeps an empty cell
    assert sorted(r["conf"]["names"]) == ["b", "c"] and r["focus"] == "an"
    assert r["both"]["names"] == ["b"] and r["sirna"]["names"] == ["c"]
    assert set(r["conf"]["chips"]) == {"rpt-srna-an conflicts"}


@pytest.mark.skipif(shutil.which("osascript") is None, reason="needs osascript")
def test_plain_list_has_no_annotation_column_or_filter(tmp_path):
    from test_html_report_js import run_js, probe_for, page
    _, text, d = page(tmp_path, [mk("a"), mk("b", reads=400)])
    r = run_js(tmp_path, text, probe_for(d, LIST_PROBE))
    assert not r["hasAn"] and r["first"]["chips"] == [] and set(r["first"]["cells"]) == {6}
    assert len(r["first"]["ths"]) == 6
    assert "Annotation" not in r["first"]["shown"]
