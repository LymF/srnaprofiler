"""A class must be earned by evidence, not by the size window alone.

Size is a filter, not a signature: at 21nt every plant locus fits phasiRNA and
siRNA equally, so whichever carries the larger w_size wins by default. That is
how 184 annotated Arabidopsis miRNAs ended up with a confidently asserted class.
With `min_evidence` set, a locus whose winning class rests only on its size is
labelled `unclassified` instead.
"""
from records import ClusterRecord
from decider import RuleScorer

CFG = {
    "decision": {"review_confidence": 0.4, "min_margin": 0.0, "min_reads": 10,
                 "min_evidence": 0.3},
    "classes": {
        "phasiRNA": {"size": {"lo": 21, "hi": 22}, "w_size": 1.0, "w_phase": 3.0,
                     "phase_x0": 5.0, "phase_scale": 1.5},
        "hcsiRNA": {"size": {"lo": 23, "hi": 24}, "w_size": 1.0, "w_5a": 1.0},
    },
}


def _decide(cfg=None, **kw):
    base = dict(chrom="c", start=0, end=100, reads=500, frac_top=0.5, pct_1u=0.1,
                pct_5a=0.1, pingpong_frac=0.0, phase_score=0.0)
    base.update(kw)
    r = ClusterRecord(**base)
    RuleScorer(cfg or CFG).decide(r)
    return r


def test_size_alone_is_not_enough():
    r = _decide(dicer_call="21")           # fits phasiRNA, no phasing at all
    assert r.cls == "unclassified"
    assert r.review is True


def test_real_evidence_earns_the_class():
    r = _decide(dicer_call="21", phase_score=30.0)
    assert r.cls == "phasiRNA"
    assert r.review is False


def test_evidence_from_another_feature_also_counts():
    r = _decide(dicer_call="24", pct_5a=0.8)
    assert r.cls == "hcsiRNA"


def test_absent_min_evidence_keeps_the_old_behaviour():
    cfg = {"decision": {k: v for k, v in CFG["decision"].items() if k != "min_evidence"},
           "classes": CFG["classes"]}
    r = _decide(cfg=cfg, dicer_call="21")
    assert r.cls == "phasiRNA"


def test_unclassified_records_what_was_missing():
    r = _decide(dicer_call="21")
    assert "evidence" in r.evidence or "size" in r.evidence
