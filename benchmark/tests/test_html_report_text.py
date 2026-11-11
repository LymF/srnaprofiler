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
    st = rt.evidence_status(mk())
    assert st["1u"] == st["10a"] == st["pingpong"] == st["phasing"] == "unavailable"


def test_explain_piRNA_names_pingpong():
    x = rt.explain(mk(e={"pingpong": 0.9, "1u": 0.9, "5a": 0.8}))
    pp = [c for c in x["checks"] if c["key"] == "pingpong"][0]
    assert pp["status"] == "present"
    assert x["verdict"]


def test_explain_ambiguous_lists_missing():
    x = rt.explain(mk(cls="ambiguous", e={"1u": 0.0, "pingpong": 0.0, "phase": 0.0, "5a": 0.0}))
    assert x["missing"]


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
    assert any("20" in x for x in w)


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
