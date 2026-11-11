import json
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import report_html
import report_text as rt
from report_data import Cluster, Selection, select

_E = ("1u", "5a", "fractop_bid", "pingpong", "phase")


def mk(name="c1", reads=500, cls="piRNA", e=None, **kw):
    ev = {k: None for k in _E}
    ev.update(e if e is not None else {k: 0.7 for k in _E})
    base = dict(name=name, chrom="chr1", start=1000, end=1100, reads=reads, reads_rpm=5.0,
                cls=cls, confidence=0.9, dicer_call=26, size_purity=0.8, frac_top=0.8,
                pct_1u=0.8, pct_10a=0.5, pingpong_z=3.0, phase_pvalue=0.5,
                multimap_frac=0.1, unique_frac=0.9, review=False, major_rna="UGAUCG",
                evidence="", e=ev, gates={"piRNA": True, "miRNA": False},
                has_evidence=any(v is not None for v in ev.values()))
    base.update(kw)
    return Cluster(**base)


def zero_sig():
    return {"sizes": {"lo": 18, "plus": [0.0] * 18, "minus": [0.0] * 18},
            "nt": [{"A": 0, "C": 0, "G": 0, "U": 0} for _ in range(26)],
            "coverage": {"start": 1000, "bin": 1, "plus": [0.0] * 100, "minus": [0.0] * 100},
            "pingpong": [0.0] * 30, "n_alns": 0}


SUMMARY = {"library_reads": 1e6, "cluster_reads": 5000.0, "unique_share": 0.8,
           "fragments": 10, "fragment_reads": 40.0}


def build(tmp_path, clusters, floor=300.0, summary=SUMMARY, title="t", max_detail=3000):
    sel = select(clusters, floor, max_detail=max_detail)
    sig = {c.name: zero_sig() for c in sel.detail}
    out = tmp_path / "r.html"
    report_html.build_report(sel, summary, sig, title, str(out))
    return sel, out.read_text(), out


def data_of(html):
    m = re.search(r'<script type="application/json" id="report-data">(.*?)</script>', html, re.S)
    return json.loads(m.group(1))


def test_output_is_self_contained(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    for bad in ('src="http', 'href="http', "<script src", "@import", "url(http"):
        assert bad not in html


def test_embedded_json_roundtrips(tmp_path):
    sel, html, _ = build(tmp_path, [mk(f"c{i}", reads=400 + i) for i in range(7)])
    d = data_of(html)
    assert len(d["rows"]) == len(sel.shown) == 7
    for k in ("library", "counts", "warnings", "rows", "details", "strings", "glossary", "references"):
        assert k in d
    assert d["counts"] == {"piRNA": 7}


def test_script_tag_cannot_be_closed(tmp_path):
    evil = '</script><b>x'
    _, html, _ = build(tmp_path, [mk(evil), mk('a&b"\'<i>')], title='</script><title>')
    assert html.count("<script") == html.count("</script>") == 2
    d = data_of(html)
    assert evil in [r["n"] for r in d["rows"]]
    assert "</script><title>" not in html
    assert "<title></script>" not in html


def test_empty_selection_builds(tmp_path):
    sel, html, _ = build(tmp_path, [mk(reads=10)], floor=300.0, summary={})
    d = data_of(html)
    assert d["rows"] == [] and d["counts"] == {}
    assert d["warnings"]
    assert "No clusters passed" in d["meta"]["summary"]


def test_old_tsv_builds(tmp_path):
    cl = [mk(f"c{i}", e={}, pct_10a=None, has_evidence=False) for i in range(3)]
    _, html, _ = build(tmp_path, cl)
    d = data_of(html)
    # statuses for 1u, 10a, pingpong, phasing are "unavailable" (digit 3)
    for r in d["rows"]:
        assert r["v"][1] == r["v"][2] == r["v"][3] == r["v"][4] == "3"


def test_size_budget(tmp_path):
    cl = [mk(f"cluster_{i:05d}", reads=300 + i) for i in range(1625)]
    _, html, out = build(tmp_path, cl)
    assert out.stat().st_size < 6_000_000


def test_cap_rows_without_detail(tmp_path):
    sel, html, _ = build(tmp_path, [mk(f"c{i}", reads=400 + i) for i in range(5)], max_detail=2)
    d = data_of(html)
    assert sum(1 for r in d["rows"] if r["d"]) == 2
    assert len(d["details"]) == 2
    assert d["meta"]["capped"] is True


def test_footer_states_display_thresholds(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    assert "display-only" in html and str(rt.DISPLAY_PRESENT) in html
