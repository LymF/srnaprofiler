"""The class score is the weakest of its defining evidence.

Minimum, not sum: summing would let one strong signal cover for a missing one,
which is the false positive being targeted. Score 2 therefore reads as "every
defining evidence has p <= 0.01".
"""
from records import ClusterRecord
from decider import RuleScorer

CFG = {
    "decision": {"review_confidence": 0.4, "min_margin": 0.0, "min_reads": 10},
    "classes": {
        "hcsiRNA": {"size": {"lo": 23, "hi": 24}, "w_size": 1.0, "w_5a": 1.0,
                    "defining_evidence": ["5a", "density"]},
        "phasiRNA": {"size": {"lo": 21, "hi": 22}, "w_size": 1.0, "w_phase": 3.0,
                     "defining_evidence": ["phase"]},
    },
}


def _rec(**kw):
    base = dict(chrom="c", start=0, end=100, name="n", reads=500,
                frac_top=0.5, pct_1u=0.1, pct_5a=0.8, pingpong_frac=0.0,
                phase_score=0.0)
    base.update(kw)
    return ClusterRecord(**base)


def test_score_is_the_minimum_of_defining_evidence():
    r = _rec(dicer_call="24", s_5a=14.0, s_density=9.0)
    RuleScorer(CFG).decide(r)
    assert r.cls == "hcsiRNA"
    assert r.score == 9.0


def test_one_missing_evidence_sinks_the_score():
    # spectacular density cannot cover for absent base bias
    r = _rec(dicer_call="24", s_5a=0.0, s_density=40.0)
    RuleScorer(CFG).decide(r)
    assert r.cls == "hcsiRNA"
    assert r.score == 0.0


def test_non_defining_evidence_is_ignored_by_the_score():
    r = _rec(dicer_call="24", s_5a=14.0, s_density=9.0, s_pingpong=99.0)
    RuleScorer(CFG).decide(r)
    assert r.score == 9.0


def test_class_without_defining_evidence_scores_zero():
    cfg = {"decision": CFG["decision"],
           "classes": {"x": {"size": {"lo": 20, "hi": 30}, "w_size": 1.0}}}
    r = _rec(dicer_call="24")
    RuleScorer(cfg).decide(r)
    assert r.score == 0.0


def test_score_belongs_to_the_winning_class():
    r = _rec(dicer_call="21", s_phase=7.0, s_5a=30.0, s_density=30.0)
    RuleScorer(CFG).decide(r)
    assert r.cls == "phasiRNA"
    assert r.score == 7.0


def test_unclassified_by_min_evidence_gate_has_zero_score():
    # winning class has defining evidence, but its support fails min_evidence,
    # so cls is downgraded to unclassified -- score must not still describe
    # the rejected candidate class.
    cfg = {
        "decision": {"review_confidence": 0.4, "min_margin": 0.0,
                     "min_reads": 10, "min_evidence": 100.0},
        "classes": CFG["classes"],
    }
    r = _rec(dicer_call="24", s_5a=14.0, s_density=9.0)
    RuleScorer(cfg).decide(r)
    assert r.cls == "unclassified"
    assert r.score == 0.0
