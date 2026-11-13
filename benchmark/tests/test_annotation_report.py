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
    assert LEGEND in data_of(html)["meta"]["an"]
    _, plain, _ = build(tmp_path, [mk("a")])
    d = data_of(plain)
    assert "an" not in d["meta"] and LEGEND not in plain
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
    _, text, d = page(tmp_path, [mk("a")])
    r = run_js(tmp_path, text, probe_for(d, probe))
    assert LEGEND not in r["detail"] and r["chips"] == []
