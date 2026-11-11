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


def test_explain_piRNA_names_pingpong():
    x = rt.explain(mk(e={"pingpong": 0.9, "1u": 0.9, "5a": 0.8}))
    pp = [c for c in x["checks"] if c["key"] == "pingpong"][0]
    assert pp["status"] == "present"
    assert x["verdict"]


def test_explain_ambiguous_gives_reason():
    x = rt.explain(mk(cls="ambiguous", e={"1u": 0.0, "pingpong": 0.0, "phase": 0.0}, pct_10a=0.0))
    assert x["unresolved"]


def test_explain_unknown_class_is_generic():
    for cls in ("other", "phasiRNA24", "hcsiRNA", "weird"):
        x = rt.explain(mk(cls=cls))
        assert x["verdict"]
        assert len(x["checks"]) == 7


def test_why_not_uses_gates():
    x = rt.explain(mk(gates={"piRNA": True, "miRNA": False}))
    assert [w["cls"] for w in x["why_not"]] == ["miRNA"]


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


def test_every_cited_key_has_reference():
    for g in rt.GLOSSARY:
        for k in g["cite"]:
            assert k in rt.REFERENCES


def test_strings_have_no_markup():
    assert all("<" not in v for v in rt.STRINGS.values())
    assert all("<" not in g["text"] for g in rt.GLOSSARY)


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

def _check(x, key):
    return [c for c in x["checks"] if c["key"] == key][0]


def test_multimapping_check_reads_high_as_bad():
    hi = _check(rt.explain(mk(multimap_frac=0.8)), "multimap")
    lo = _check(rt.explain(mk(multimap_frac=0.1)), "multimap")
    mid = _check(rt.explain(mk(multimap_frac=rt.DISPLAY_MULTIMAP_MODERATE)), "multimap")
    assert (hi["status"], lo["status"], mid["status"]) == ("high", "low", "moderate")
    assert rt.STRINGS["check.multimap.title"] == "Multimapping"
    assert "Many reads map" in " ".join(rt.explain(mk(multimap_frac=rt.DISPLAY_MULTIMAP_HIGH))["cautions"])
    assert hi["coloured"]


def test_low_confidence_caution_follows_review_flag():
    assert not rt.explain(mk(confidence=0.2, has_evidence=True))["cautions"]
    assert rt.STRINGS["caution.review"] in rt.explain(mk(review=True))["cautions"]
    assert "caution.lowconf" not in rt.STRINGS


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
    assert crit["unresolved"] in crit["verdict"]
    assert rt.explain(mk())["unresolved"] == ""


def test_ambiguous_text_does_not_say_fits_several_classes():
    amb = [g for g in rt.GLOSSARY if g["term"] == "ambiguous"][0]["text"]
    for t in (amb, rt.STRINGS["verdict.ambiguous"], rt.STRINGS["ui.stat.amb.sub"],
              rt.STRINGS["ui.help.classes"]):
        assert "more than one class" not in t and "single class" not in t
    assert "review" in amb


def test_confidence_text_states_floor():
    conf = [g for g in rt.GLOSSARY if g["term"] == "confidence"][0]["text"]
    assert "from 0 to 1" not in conf and "1 divided by" in conf
    assert "{min}" in rt.STRINGS["ui.detail.conf"]
    assert rt.confidence_floor("animal") == 0.5
    assert rt.confidence_floor("plant") == 0.25


def test_defining_checks_coloured_others_context():
    x = rt.explain(mk(cls="piRNA"))
    col = {c["key"] for c in x["checks"] if c["coloured"]}
    assert col == {"size", "1u", "pingpong", "multimap"}
    s = rt.explain(mk(cls="siRNA"))
    assert {c["key"] for c in s["checks"] if c["coloured"]} == {"size", "strand", "multimap"}
    a = rt.explain(mk(cls="ambiguous", evidence="nearest=piRNA; x"))
    assert {c["key"] for c in a["checks"] if c["coloured"]} == {"multimap"}
    assert a["note"]
    m = rt.explain(mk(cls="miRNA"))
    assert "miRDeep2" in m["note"]
    h = rt.explain(mk(cls="hcsiRNA"), kingdom="plant")
    assert {c["key"] for c in h["checks"] if c["coloured"]} == {"size", "5a", "multimap"}
    assert not rt.explain(mk(cls="piRNA"))["note"]


def test_plant_checks_and_text_have_no_piRNA():
    x = rt.explain(mk(cls="hcsiRNA", e={"5a": 0.8}), kingdom="plant")
    keys = [c["key"] for c in x["checks"]]
    assert "5a" in keys and "10a" not in keys and "pingpong" not in keys
    assert _check(x, "5a")["status"] == "present"
    texts = [c["detail"] for c in x["checks"]] + [x["verdict"]]
    texts += [rt.STRINGS[f"verdict.{k}"] for k in ("phasiRNA", "phasiRNA24", "hcsiRNA", "siRNA", "miRNA")]
    texts += [rt.size_help("plant")] + [rt.dot_legend("plant")]
    assert not any("piRNA" in t.replace("phasiRNA", "") for t in texts)
    assert "piRNA" in rt.size_help("animal")
    assert list(rt.evidence_status(mk(), "plant")) == ["size", "1u", "5a", "strand", "phasing"]
    assert list(rt.evidence_status(mk())) == ["size", "1u", "10a", "pingpong", "phasing"]


def test_hcsiRNA_glossary_cites_verified_sources():
    g = [g for g in rt.GLOSSARY if g["term"] == "hcsiRNA"][0]
    assert "mi2008" in g["cite"] and "10.1016/j.cell.2008.02.034" in rt.REFERENCES["mi2008"]


def test_why_not_names_failing_criterion():
    x = rt.explain(mk(cls="piRNA", dicer_call=26, gates={"piRNA": True, "siRNA": False}))
    w = x["why_not"][0]
    assert w["cls"] == "siRNA" and "26" in w["reason"] and "20" in w["reason"] and "22" in w["reason"]
    p = rt.explain(mk(cls="hcsiRNA", dicer_call=24, phase_pvalue=0.4,
                      gates={"phasiRNA24": False, "hcsiRNA": True}), kingdom="plant")
    assert "0.4" in p["why_not"][0]["reason"]
    short = rt.explain(mk(cls="siRNA", dicer_call=26, start=0, end=300,
                          gates={"piRNA": False, "siRNA": True}))
    assert "1000" in short["why_not"][0]["reason"]


def test_why_not_fallback_when_nothing_failed():
    x = rt.explain(mk(gates={"piRNA": True, "siRNA": True}))
    assert x["why_not"] == [{"cls": "", "reason": rt.STRINGS["whynot.none"]}]


def test_minor_wording():
    assert "within 1 nt" in rt.STRINGS["check.size.detail"]
    assert "not weighted" in rt.STRINGS["ui.fig.cov.help"]
    assert "--min-reads" in rt.STRINGS["warn.few"]
