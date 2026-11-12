import json
import os
import re
import sys

import pytest

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
    z = lambda n: [0.0] * n
    return {"n_alns": 0, "dominant": None,
            "sizes": dict({"lo": 18}, **{b: {"+": z(18), "-": z(18)} for b in "ACGU"}),
            "phasing": {"+": z(100), "-": z(100)}, "pingpong": z(30),
            "cov": {"start": 1000, "bin": 1, "+": z(100), "-": z(100)},
            "ends": {"bin": 1, "+": z(100), "-": z(100)},
            "logo": {s: [dict.fromkeys("ACGU", 0)] * 26 for s in "+-"}, "seqs": []}


def full_sig(rng=None, bins=400):
    """Task-1-shaped signal with non-zero random magnitudes and full lengths."""
    import random
    r = rng or random.Random(1)
    arr = lambda n, hi: [round(r.random() * hi * (r.random() < .7), 3) for _ in range(n)]
    lg = [dict(zip("ACGU", [round(x / sum(w), 3) for x in w]))
          for w in ([r.random() + .1 for _ in range(4)] for _ in range(26))]
    seq = lambda: "".join(r.choice("ACGU") for _ in range(r.randint(21, 31)))
    sq = [[q, len(q), r.choice("+-"), round(r.random() * 5000, 3), q[0]]
          for q in (seq() for _ in range(10))]
    return {"n_alns": 900, "dominant": 24,
            "sizes": dict({"lo": 18}, **{b: {s: arr(18, 800) for s in "+-"} for b in "ACGU"}),
            "phasing": {s: arr(100, 300) for s in "+-"}, "pingpong": arr(30, 500),
            "cov": {"start": 1000, "bin": 3, "+": arr(bins, 90), "-": arr(bins, 90)},
            "ends": {"bin": 3, "+": arr(bins, 60), "-": arr(bins, 60)},
            "logo": {"+": lg, "-": lg}, "seqs": sq}


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
    for k in ("library", "counts", "warnings", "rows", "details", "strings"):
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


def test_cap_rows_without_detail(tmp_path):
    sel, html, _ = build(tmp_path, [mk(f"c{i}", reads=400 + i) for i in range(5)], max_detail=2)
    d = data_of(html)
    assert sum(1 for r in d["rows"] if r["d"]) == 2
    assert len(d["details"]) == 2
    assert d["meta"]["capped"] is True


def test_footer_states_display_thresholds(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    assert "display-only" in html and str(rt.DISPLAY_PRESENT) in html


_APP = open(os.path.join(os.path.dirname(__file__), "..", "..", "report_assets", "app.js"),
            encoding="utf-8").read()
_CODE = re.sub(r"(?m)^\s*//.*$", "", _APP)


@pytest.mark.xfail(strict=False, reason="front end lands in Task 5")
def test_app_text_keys_exist_in_strings():
    keys = re.findall(r'\bT\("([^"]+)"', _CODE)
    assert keys
    missing = [k for k in keys if not k.endswith(".") and k not in rt.STRINGS]
    assert not missing
    # dynamic keys: check.<key>.title and status.<level>
    for k in ("size", "1u", "5a", "10a", "strand", "pingpong", "phasing", "multimap"):
        assert f"check.{k}.title" in rt.STRINGS
    assert "check.context" in rt.STRINGS
    for lv in ("present", "weak", "absent", "unavailable", "low", "moderate", "high"):
        assert f"status.{lv}" in rt.STRINGS
    for k in re.findall(r'"(ui\.[\w.]+)"', _CODE):
        assert k in rt.STRINGS


def test_app_has_no_hardcoded_english_phrases():
    # A literal is a phrase when it is only letters, spaces and light punctuation and has
    # two or more words (CSS class lists contain '-', '_' or digits and are skipped).
    # Allowed: "use strict" and the CSS pairs "dot <status>" / "ico <status>".
    pat = re.compile(r'"([A-Za-z][A-Za-z\' ,.:;!?]*[A-Za-z]) ?"')
    bad = []
    for m in pat.finditer(_CODE):
        lit = m.group(1).strip()
        if " " not in lit or lit == "use strict" or lit.split()[0] in ("dot", "ico"):
            continue
        bad.append(lit)
    assert bad == []


def test_footer_and_title_come_from_strings(tmp_path):
    _, html, _ = build(tmp_path, [mk()], title="")
    assert "display-only" in rt.STRINGS["ui.footer"]
    assert rt.STRINGS["ui.title.default"] in html


# ---- final review fixes ----

def test_ambiguous_row_has_no_confidence_and_comes_last(tmp_path):
    amb = mk("amb", reads=900, cls="ambiguous", confidence=0.99,
             evidence="nearest=piRNA; size=0.6", dicer_call=23)
    _, html, _ = build(tmp_path, [amb, mk("ok", confidence=0.6)])
    d = data_of(html)
    assert [r["n"] for r in d["rows"]] == ["ok", "amb"]
    r = d["rows"][1]
    assert r["f"] is None and r["nc"] == "piRNA" and "23" in r["u"]
    assert "nc" not in d["rows"][0] and "u" not in d["rows"][0]


def test_strip_levels_and_context(tmp_path):
    _, html, _ = build(tmp_path, [mk(multimap_frac=0.8)])
    det = data_of(html)["details"]["0"]
    st = {k: (lv, v) for k, lv, v in det["strip"]}
    assert st["multimap"] == (report_html._LEVEL["high"], 0.8)
    assert st["size"][1] == 26 and st["pingpong"][1] == 3.0 and st["1u"][1] == 0.8
    assert det["v"] == "verdict.piRNA" and det["v"] in rt.STRINGS
    assert det["ctx"] == [100, None]


def test_strip_has_measured_values_for_old_tsv(tmp_path):
    cl = mk(e={}, pct_1u=None, pct_10a=None, pingpong_z=None, phase_pvalue=None,
            frac_top=None, multimap_frac=None, dicer_call=None, has_evidence=False)
    _, html, _ = build(tmp_path, [cl])
    st = {k: (lv, v) for k, lv, v in data_of(html)["details"]["0"]["strip"]}
    for k in ("1u", "10a", "pingpong", "phasing", "strand", "multimap"):
        assert st[k] == (report_html._LEVEL["unavailable"], None), k


def test_payload_has_no_glossary_or_references(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    d = data_of(html)
    assert "glossary" not in d and "references" not in d
    assert not any(k.startswith(("ui.gl.", "ui.ref.")) for k in d["strings"])
    for k in ("checks", "note", "why", "cautions"):
        assert k not in d["details"]["0"]


def _worst(pair, orig):
    m, q = pair
    got = report_html.dequantise(pair, len(orig))
    return max(abs(a - b) for a, b in zip(got, orig)), m


def test_quantise_roundtrip_error_bound():
    import random
    r = random.Random(3)
    for _ in range(50):
        orig = [r.random() * 10 ** r.randint(-2, 4) * (r.random() < .8) for _ in range(r.randint(1, 400))]
        pair = report_html._quant(orig)
        err, m = _worst(pair, orig)
        assert err <= m / 198 * 1.001 + 1e-9
        assert all(0 <= i <= 99 for i in pair[1]) and (not pair[1] or pair[1][-1] != 0)


def test_quantise_all_zero_array():
    assert report_html._quant([0.0] * 30) == [0, []]
    assert report_html._quant([]) == [0, []]
    assert report_html._quant([float("nan"), 0]) == [0, []]
    assert report_html.dequantise([0, []], 5) == [0.0] * 5


def test_signals_quantised_shape_and_bound(tmp_path):
    sig = full_sig()
    out = report_html._signals(sig)
    m, q = out["cv"]["+"]
    assert len(q) <= 400 and out["sq"] == sig["seqs"] and out["dom"] == 24
    err, mx = _worst(out["ph"]["+"], sig["phasing"]["+"])
    assert err <= mx / 198 * 1.001
    assert set(out["sz"]) == {"lo", "A", "C", "G", "U"}
    assert all(0 <= x <= 99 for pos in out["lg"]["+"] for x in pos) and len(out["lg"]["+"]) == 26
    json.dumps(out, allow_nan=False)


def test_size_budget_realistic(tmp_path):
    import random
    rng = random.Random(7)
    cl = [mk(f"cluster_{i:05d}", reads=300 + i) for i in range(1000)]
    sel = select(cl, 300.0)
    sig = {c.name: full_sig(rng) for c in sel.detail}
    out = tmp_path / "r.html"
    report_html.build_report(sel, SUMMARY, sig, "t", str(out))
    d = data_of(out.read_text())
    per = [len(json.dumps(v, separators=(",", ":"))) for v in d["details"].values()]
    print("max per-cluster detail JSON bytes", max(per), "file", out.stat().st_size)
    assert max(per) < 9_000
    assert out.stat().st_size < 12_000_000


def test_about_strings_in_table(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    d = data_of(html)
    S = d["strings"]
    for k in ("ui.about.title", "ui.about.tagline", "ui.about.github",
              "ui.about.install.1", "ui.about.install.2"):
        assert S[k] and "<" not in S[k]
    m = d["meta"]
    import version
    assert m["version"] == version.__version__ and m["github"].startswith("https://github.com/")
    assert m["install"] == [S["ui.about.install.1"], S["ui.about.install.2"]]
    assert m["tagline"] == S["ui.about.tagline"]


def test_name_with_script_tag_still_escaped(tmp_path):
    _, html, _ = build(tmp_path, [mk("</script><b>x")])
    assert html.count("</script>") == 2 and "</script><b>x" not in html


def test_meta_is_kingdom_aware(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    m = data_of(html)["meta"]
    assert m["kingdom"] == "animal" and m["conf_min"] == 0.5 and m["conf_k"] == 2
    assert "piRNA" in m["kt"]["sizes"] and m["dot_names"][2] == rt.STRINGS["ui.dot.10a"]
    p = mk(cls="hcsiRNA", gates={"hcsiRNA": True, "siRNA": False})
    _, html, _ = build(tmp_path, [p])
    d = data_of(html)
    assert d["meta"]["kingdom"] == "plant" and d["meta"]["conf_min"] == 0.25
    assert "piRNA" not in d["meta"]["kt"]["sizes"].replace("phasiRNA", "")
    assert [c[0] for c in d["details"]["0"]["strip"]][2] == "5a"
    assert d["meta"]["defining"]["hcsiRNA"] == ["5a", "size"]


def test_no_reads_and_failed_markers(tmp_path):
    sel = select([mk("a"), mk("b"), mk("c")], 0)
    full = zero_sig()
    full["n_alns"] = 3
    sig = {"a": zero_sig(), "b": dict(zero_sig(), failed=True), "c": full}
    out = tmp_path / "r.html"
    report_html.build_report(sel, SUMMARY, sig, "t", str(out))
    det = data_of(out.read_text())["details"]
    by = {sel.shown[int(i)].name: v["sig"] for i, v in det.items()}
    assert by["a"] == {"none": 1} and by["b"] == {"err": 1}
    assert "cv" in by["c"] and "sz" in by["c"]


def test_footer_lists_multimap_cutoffs(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    assert "Multimapping: moderate at 0.3" in html and "high at 0.6" in html


def test_meta_floor_is_json_safe(tmp_path):
    sel = select([mk()], 0)
    sel.floor = float("nan")
    out = tmp_path / "r.html"
    report_html.build_report(sel, SUMMARY, {}, "t", str(out))
    assert data_of(out.read_text())["meta"]["floor"] is None


@pytest.mark.xfail(strict=False, reason="front end lands in Task 5")
def test_mirna_confidence_note():
    # classify sets confidence 1.0 for miRNA (miRDeep2), not a softmax over classes
    assert 'T("ui.detail.conf.mirna")' in _CODE
    assert "miRDeep2" in rt.STRINGS["ui.detail.conf.mirna"]
