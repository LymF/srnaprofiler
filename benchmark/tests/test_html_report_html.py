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
        if (" " not in lit or lit in ("use strict", "noopener noreferrer")
                or lit.split()[0] in ("dot", "ico")):
            continue
        bad.append(lit)
    assert bad == []


def test_footer_and_title_come_from_strings(tmp_path):
    _, html, _ = build(tmp_path, [mk()], title="")
    assert all("display-only" in rt.STRINGS[f"ui.footer.{k}"] for k in ("animal", "plant"))
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
    assert set(out["sz"]) == {"lo", "max", "A", "C", "G", "U"}
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


def test_mirna_confidence_note():
    # classify sets confidence 1.0 for miRNA (miRDeep2), not a softmax over classes
    assert 'T("ui.detail.conf.mirna")' in _CODE
    assert "miRDeep2" in rt.STRINGS["ui.detail.conf.mirna"]


def test_signals_carry_read_counts_and_bins():
    sig = full_sig()
    out = report_html._signals(sig)
    assert out["n"] == 900 and out["nb"] == 400
    want = sum(sig["sizes"][b][s][24 - 18] for b in "ACGU" for s in "+-")
    assert abs(out["nd"] - want) < 0.06
    assert report_html._signals(dict(sig, dominant=None))["nd"] == 0


def test_sizes_stack_error_within_one_percent_of_panel_max():
    import random
    r = random.Random(11)
    for _ in range(30):
        sig = full_sig(r)
        for b in "ACGU":
            for s in "+-":
                sig["sizes"][b][s] = [v * r.random() for v in sig["sizes"][b][s]]
        sz = report_html._signals(sig)["sz"]
        top = max(max(sig["sizes"][b][s]) for b in "ACGU" for s in "+-")
        assert abs(sz["max"] - top) <= top * 1e-5
        for si, s in enumerate("+-"):
            for i in range(18):
                got = sum(report_html.dequantise([sz["max"], sz[b][si]], 18, top=999)[i]
                          for b in "ACGU")
                want = sum(sig["sizes"][b][s][i] for b in "ACGU")
                assert abs(got - want) <= 0.01 * sz["max"]
    assert report_html._quant_group([[0.0], [0.0]]) == (0, [[], []])


def test_quant_keeps_scale_of_tiny_max():
    m, q = report_html._quant([1e-12, 0])
    assert m == 1e-12 and q == [99]


def test_meta_defining_animal(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    df = data_of(html)["meta"]["defining"]
    assert df["piRNA"] == ["1u", "pingpong", "size"]
    assert df["siRNA"] == ["size", "strand"]
    assert "hcsiRNA" not in df


def test_logo_tolerates_nan_and_none():
    sig = full_sig()
    sig["logo"]["+"][0] = {"A": float("nan"), "C": None, "G": 0.5, "U": 0.5}
    out = report_html._signals(sig)
    assert out["lg"]["+"][0] == [0, 0, 50, 50]
    json.dumps(out, allow_nan=False)


def test_review_flag_reaches_rows(tmp_path):
    _, html, _ = build(tmp_path, [mk("a", reads=900, review=True), mk("b", reads=500)])
    rows = {r["n"]: r for r in data_of(html)["rows"]}
    assert rows["a"]["rv"] == 1 and "rv" not in rows["b"]


# ---- v2 front end: fonts, motion, text size ----

_ASSETS = os.path.join(os.path.dirname(__file__), "..", "..", "report_assets")
_CSS = open(os.path.join(_ASSETS, "report.css"), encoding="utf-8").read()
_MARK = "srnaprofiler additions"
_OWN_CSS = _CSS[_CSS.index(_MARK):]


def test_fonts_are_inlined_and_no_inter(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    faces = re.findall(r"@font-face\{([^}]*)\}", html)
    assert 1 <= len(faces) <= 4
    fams = {re.search(r"font-family:'([^']+)'", f).group(1) for f in faces}
    assert fams == {"IBM Plex Sans", "IBM Plex Mono"}
    assert {re.search(r"font-weight:(\d+)", f).group(1) for f in faces} == {"400", "500"}
    assert all("url(data:font/woff2;base64," in f for f in faces)
    assert "fonts.googleapis" not in html and "fonts.gstatic" not in html
    assert not re.search(r"\bInter\b", html)
    assert len([f for f in os.listdir(os.path.join(_ASSETS, "fonts")) if f.endswith(".woff2")]) <= 4
    assert os.path.exists(os.path.join(_ASSETS, "fonts", "LICENSE-IBM-Plex.txt"))
    # weights used by the stylesheet's own block stay at 400/500
    assert set(re.findall(r"font-weight:\s*(\d+)", _OWN_CSS)) <= {"400", "500"}


def test_no_scroll_reveal():
    for src in (_APP, _CSS):
        for bad in ("IntersectionObserver", "animation-timeline", "scroll-timeline",
                    '"scroll"', "'scroll'", "onscroll"):
            assert bad not in src, bad
        banned = {"all", "width", "height", "top", "left", "right", "bottom", "margin", "padding",
                  "max-height", "min-height", "inset"}
        for decl in re.findall(r"transition(?:-property)?\s*:([^;]*);", src):
            props = {p.split()[0] for p in decl.split(",") if p.split()}
            assert not props & banned and not any(p.startswith(("margin-", "padding-"))
                                                  for p in props), decl
    # keyframes of the own block move only opacity and transform
    for body in re.findall(r"@keyframes[^{]*\{(.*?\})\s*\}", _OWN_CSS, re.S):
        props = set(re.findall(r"([a-z-]+)\s*:", body))
        assert props <= {"opacity", "transform"}, props


def test_reduced_motion_block_exists():
    m = re.search(r"@media \(prefers-reduced-motion: reduce\)\s*\{(.*?)\n\}", _OWN_CSS, re.S)
    assert m
    body = m.group(1)
    assert "animation: none" in body and "transition: none" in body
    assert "opacity: 1" in body                      # nothing stays hidden
    for tok in ("--m-in: 150ms", "--m-out: 120ms", "--m-press: 90ms", "--m-panel: 200ms",
                "--m-slide: 240ms", "cubic-bezier(.2, .8, .2, 1)", "cubic-bezier(.22, 1, .36, 1)"):
        assert tok in _OWN_CSS, tok


def test_min_text_size():
    sizes = [float(x) for x in re.findall(r"font-size:\s*([\d.]+)px", _OWN_CSS)]
    assert sizes and min(sizes) >= 13
    rems = [float(x) for x in re.findall(r"font-size:\s*([\d.]+)rem", _OWN_CSS)]
    assert all(r * 16 >= 13 for r in rems)
    for tok in ("--rpt-text-xs", "--rpt-text-sm"):
        v = float(re.search(tok + r":\s*([\d.]+)rem", _CSS).group(1))
        assert v * 16 >= 13
    # the SVG helpers set one size constant, never a literal under 13
    fs = re.search(r"\bvar FS = (\d+);", _CODE)
    assert fs and int(fs.group(1)) >= 13
    for lit in re.findall(r'"font-size"\s*:\s*([\d.]+)', _CODE):
        assert float(lit) >= 13
    assert not re.search(r"font-?[sS]ize[^;\n]*\b([0-9]|1[0-2])px", _CODE)


def test_app_js_has_no_glossary_code():
    for bad in ("glossary", "Glossary", "references", "renderRefs", "refNum", "cite",
                "ui.gl.", "ui.ref.", "ui.howto", "ui.detail.why", "whynot", "cautions",
                "d.checks", "d.why"):
        assert bad not in _CODE, bad


def test_app_has_no_inner_html():
    for bad in ("innerHTML", "outerHTML", "insertAdjacentHTML", "document.write"):
        assert bad not in _APP


def test_strings_have_no_orphans():
    code = _CODE
    py = ""
    for f in ("report_text.py", "report_html.py", "report_mirna.py"):
        src = open(os.path.join(os.path.dirname(__file__), "..", "..", f), encoding="utf-8").read()
        # the STRINGS table itself does not count as a use
        py += re.sub(r"STRINGS = \{.*?\n\}", "", src, flags=re.S)
    dynamic = ("verdict.", "status.", "check.", "amb.crit.", "ui.dot.", "ui.help.sizes.",
                   "ck.thr.", "ck.obs.", "ck.st.", "lq.", "mi.nostruct.", "mi.num.", "mi.crit.")
    orphans = [k for k in rt.STRINGS
               if f'"{k}"' not in code and f'"{k}"' not in py and not k.startswith(dynamic)]
    assert orphans == []


def _call_vars(code, start):
    """Keys of the object literal passed after T("key", ...), or None when there is none."""
    i = start
    while code[i] in " \n":
        i += 1
    if code[i] != ",":
        return None
    i = code.index("{", i)
    depth, j = 0, i
    while True:
        depth += {"{": 1, "}": -1}.get(code[j], 0)
        if depth == 0:
            break
        j += 1
    body, flat, depth = code[i + 1:j], "", 0
    for ch in body:                                # keep only the top level of the literal
        depth += {"{": 1, "(": 1, "[": 1, "}": -1, ")": -1, "]": -1}.get(ch, 0)
        flat += ch if depth == 0 or ch in "{([" and depth == 1 else " "
    return set(re.findall(r"(\w+)\s*:", flat))


def test_placeholders_match_their_calls():
    for m in re.finditer(r'\bT\("([^"]+)"(?!\s*\+)', _CODE):
        key = m.group(1)
        need = set(re.findall(r"\{(\w+)\}", rt.STRINGS[key]))
        got = _call_vars(_CODE, m.end())
        if need:
            assert got is not None and need <= got, (key, need, got)
        else:
            assert not got, (key, got)


def test_footer_is_worded_by_kingdom(tmp_path):
    _, html, _ = build(tmp_path, [mk()])
    foot = re.search(r'<footer class="rpt-srna-footer">(.*?)</footer>', html, re.S).group(1)
    assert "ping-pong" in foot and "position 10" in foot and "5&#x27;A" not in foot
    p = mk(cls="hcsiRNA", gates={"hcsiRNA": True, "siRNA": False})
    _, html, _ = build(tmp_path, [p])
    foot = re.search(r'<footer class="rpt-srna-footer">(.*?)</footer>', html, re.S).group(1)
    assert "ping-pong" not in foot and "position 10" not in foot
    for want in ("Size consistency", "1U", "5&#x27;A", "strand", "phasing", "Multimapping"):
        assert want in foot, want


def test_percentile_string_has_no_ordinal():
    assert not re.search(r"\d?\}(st|nd|rd|th)\b", rt.STRINGS["ui.ctx.pct"])
