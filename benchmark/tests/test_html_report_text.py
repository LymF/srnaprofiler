import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import report_text as rt
from report_data import Cluster, Selection

_E = ("1u", "5a", "fractop_bid", "pingpong", "phase")


def mk(cls="piRNA", purity=0.8, e=None, gates=None, **kw):
    ev = {k: None for k in _E}
    ev.update(e or {})
    base = dict(name="c1", chrom="chr1", start=1, end=100, reads=500, reads_rpm=5.0,
                cls=cls, confidence=0.9, dicer_call=0, size_purity=purity, frac_top=0.8,
                pct_1u=0.8, pct_10a=0.5, pingpong_z=3.0, phase_pvalue=0.5,
                multimap_frac=0.1, unique_frac=0.9, review=False, major_rna="",
                evidence="", e=ev, gates=gates or {}, has_evidence=bool(e))
    base.update(kw)
    return Cluster(**base)


def mksel(n, floor=5.0, capped=False):
    cl = [mk() for _ in range(n)]
    return Selection(shown=cl, detail=cl, n_total=n + 5, n_hidden=5, floor=floor, capped=capped)


def test_status_thresholds():
    assert rt.evidence_status(mk(purity=0.5))["size"] == "present"
    assert rt.evidence_status(mk(purity=0.3))["size"] == "weak"
    assert rt.evidence_status(mk(purity=0.29))["size"] == "absent"
    assert rt.evidence_status(mk(e={"1u": 0.25}))["1u"] == "weak"
    assert rt.evidence_status(mk(e={"1u": 0.5}))["1u"] == "present"
    assert rt.evidence_status(mk(e={"1u": 0.1}))["1u"] == "absent"


def test_status_unavailable_on_old_tsv():
    st = rt.evidence_status(mk(pct_10a=None))
    assert st["1u"] == st["10a"] == st["pingpong"] == st["phasing"] == "unavailable"


def _strip(x, key):
    return [r for r in x["strip"] if r[0] == key][0]


def test_explain_piRNA_strip_and_verdict():
    x = rt.explain(mk(e={"pingpong": 0.9, "1u": 0.9, "5a": 0.8}))
    assert _strip(x, "pingpong")[1] == "present" and _strip(x, "pingpong")[2] == 3.0
    assert x["verdict"] == "verdict.piRNA" and x["verdict"] in rt.STRINGS


def test_explain_ambiguous_gives_reason():
    x = rt.explain(mk(cls="ambiguous", e={"1u": 0.0, "pingpong": 0.0, "phase": 0.0}, pct_10a=0.0))
    assert x["unresolved"]


def test_explain_unknown_class_is_generic():
    for cls in ("other", "phasiRNA24", "hcsiRNA", "weird"):
        x = rt.explain(mk(cls=cls))
        assert x["verdict"] in rt.STRINGS
        assert len(x["strip"]) == 7


def test_summary_sentence_names_top_class():
    s = rt.summary_sentence({"piRNA": 10, "siRNA": 3}, 5.0)
    assert "piRNA" in s and "10" in s


def test_warning_few_clusters():
    w = rt.warnings(mksel(5), {"unique_share": 0.9})
    assert any(str(rt._MIN_CLUSTERS) in x for x in w)


def test_warning_cap():
    assert len(rt.warnings(mksel(30, capped=True), {"unique_share": 0.9})) == 1


def test_warning_low_unique():
    assert rt.warnings(mksel(30), {"unique_share": 0.1})
    assert not rt.warnings(mksel(30), {"unique_share": None})


def test_strings_have_no_markup():
    assert all("<" not in v for v in rt.STRINGS.values())
    assert not hasattr(rt, "GLOSSARY") and not hasattr(rt, "REFERENCES")


def test_10a_status_from_pct_10a():
    f = lambda v, **k: rt.evidence_status(mk(pct_10a=v, **k))["10a"]
    assert f(0.5) == "present"
    assert f(0.25) == "weak"
    assert f(0.1) == "absent"
    assert f(None) == "unavailable"
    assert f(0.1, e={"5a": 0.9}) == "absent"


def test_summary_sentence_empty_and_zero():
    assert rt.summary_sentence({}, 5.0)
    assert rt.summary_sentence({"piRNA": 0}, 5.0) == rt.summary_sentence({}, 5.0)


def test_warnings_empty_selection():
    sel = Selection(shown=[], detail=[], n_total=3, n_hidden=3, floor=1000.0, capped=False)
    w = rt.warnings(sel, None)
    assert len(w) == 1 and "0" in w[0]


# ---- final review fixes ----

def test_multimapping_strip_reads_high_as_bad():
    lv = lambda v: _strip(rt.explain(mk(multimap_frac=v)), "multimap")[1]
    assert (lv(0.8), lv(0.1), lv(rt.DISPLAY_MULTIMAP_MODERATE), lv(None)) == (
        "high", "low", "moderate", "unavailable")
    assert rt.STRINGS["check.multimap.title"] == "Multimapping"


def test_nearest_class_parsed_from_evidence():
    assert rt.nearest_class(mk(cls="ambiguous", evidence="nearest=piRNA; size=0.6")) == "piRNA"
    assert rt.nearest_class(mk(cls="ambiguous", evidence="ml:hgb")) is None
    assert rt.nearest_class(mk(cls="ambiguous", evidence="nearest=siRNA; x", dicer_call=None)) is None


def test_ambiguous_reason_from_evidence():
    crit = rt.explain(mk(cls="ambiguous", dicer_call=27,
                         evidence="nearest=piRNA; below_cluster_criterion(length,density); size=1.0"))
    assert "piRNA" in crit["unresolved"] and "length" in crit["unresolved"]
    assert rt.STRINGS["amb.crit.density"] in crit["unresolved"]
    size = rt.explain(mk(cls="ambiguous", dicer_call=19, evidence="nearest=siRNA; size=0.6"))
    assert "19" in size["unresolved"] and "siRNA" in size["unresolved"]
    nosize = rt.explain(mk(cls="ambiguous", dicer_call=None, evidence="nearest=siRNA; size=0.0"))
    assert nosize["unresolved"] == rt.STRINGS["amb.nosize"]
    ml = rt.explain(mk(cls="ambiguous", evidence="ml:hgb"))
    assert ml["unresolved"] == rt.STRINGS["amb.ml"]
    other = rt.explain(mk(cls="ambiguous", evidence="something new"))
    assert other["unresolved"] == rt.STRINGS["amb.generic"]
    assert crit["verdict"] == "verdict.ambiguous"
    assert rt.explain(mk())["unresolved"] == ""


def test_ambiguous_text_does_not_say_fits_several_classes():
    for t in (rt.STRINGS["verdict.ambiguous"], rt.STRINGS["ui.stat.amb.sub"],
              rt.STRINGS["ui.help.classes"]):
        assert "more than one class" not in t and "single class" not in t


def test_confidence_floor_and_text():
    assert "{min}" in rt.STRINGS["ui.detail.conf"]
    assert rt.confidence_floor("animal") == 0.5
    assert rt.confidence_floor("plant") == 0.25


def test_plant_strip_and_text_have_no_piRNA():
    x = rt.explain(mk(cls="hcsiRNA", e={"5a": 0.8}, pct_5a=0.7), kingdom="plant")
    keys = [r[0] for r in x["strip"]]
    assert "5a" in keys and "10a" not in keys and "pingpong" not in keys
    assert _strip(x, "5a")[1] == "present" and _strip(x, "5a")[2] == 0.7
    texts = [rt.STRINGS[f"verdict.{k}"] for k in ("phasiRNA", "phasiRNA24", "hcsiRNA", "siRNA", "miRNA")]
    texts += [rt.size_help("plant")] + [rt.dot_legend("plant")]
    assert not any("piRNA" in t.replace("phasiRNA", "") for t in texts)
    assert "piRNA" in rt.size_help("animal")
    assert list(rt.evidence_status(mk(), "plant")) == ["size", "1u", "5a", "strand", "phasing"]
    assert list(rt.evidence_status(mk())) == ["size", "1u", "10a", "pingpong", "phasing"]


def test_minor_wording():
    assert "not weighted" in rt.STRINGS["ui.fig.cov.help"]
    assert "--min-reads" in rt.STRINGS["warn.few"]
