"""The class score is the weakest of its defining evidence.

Minimum, not sum: summing would let one strong signal cover for a missing one,
which is the false positive being targeted. `score` is built from the `e_*`
effect strengths, each in [0,1], so it reads as "the weakest defining signal's
strength" -- not a p-value.
"""
import pytest

from records import ClusterRecord
from decider import RuleScorer

CFG = {
    "decision": {"review_confidence": 0.4, "min_margin": 0.0, "min_reads": 10},
    "classes": {
        "hcsiRNA": {"size": {"lo": 23, "hi": 24}, "w_size": 1.0, "w_5a": 1.0,
                    "defining_evidence": ["5a", "fractop_bid"]},
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
    r = _rec(dicer_call="24", e_5a=0.9, e_fractop_bid=0.6)
    RuleScorer(CFG).decide(r)
    assert r.cls == "hcsiRNA"
    assert r.score == 0.6


def test_one_missing_evidence_sinks_the_score():
    # spectacular strand bias cannot cover for absent base bias
    r = _rec(dicer_call="24", e_5a=0.0, e_fractop_bid=0.9)
    RuleScorer(CFG).decide(r)
    assert r.cls == "hcsiRNA"
    assert r.score == 0.0


def test_non_defining_evidence_is_ignored_by_the_score():
    r = _rec(dicer_call="24", e_5a=0.8, e_fractop_bid=0.5, e_pingpong=0.99)
    RuleScorer(CFG).decide(r)
    assert r.score == 0.5


def test_class_without_defining_evidence_scores_zero():
    cfg = {"decision": CFG["decision"],
           "classes": {"x": {"size": {"lo": 20, "hi": 30}, "w_size": 1.0}}}
    r = _rec(dicer_call="24")
    RuleScorer(cfg).decide(r)
    assert r.score == 0.0


def test_score_belongs_to_the_winning_class():
    r = _rec(dicer_call="21", e_phase=0.65, e_5a=0.9, e_fractop_bid=0.9)
    RuleScorer(CFG).decide(r)
    assert r.cls == "phasiRNA"
    assert r.score == 0.65


def test_score_is_built_from_the_strengths_not_the_p_values():
    r = _rec(dicer_call="24", s_5a=99.0, s_fractop=99.0, e_5a=0.7, e_fractop_bid=0.4)
    RuleScorer(CFG).decide(r)
    assert r.score == 0.4


def test_an_unknown_evidence_name_fails_loudly():
    bad = {"decision": CFG["decision"],
           "classes": {"hcsiRNA": {"size": {"lo": 23, "hi": 24}, "w_size": 1.0,
                                   "defining_evidence": ["5a", "fractopp"]}}}
    with pytest.raises(KeyError) as e:
        RuleScorer(bad)
    assert "fractopp" in str(e.value)


def test_a_record_without_the_strength_attribute_still_scores_zero():
    """extras/simulate_rules.py rebuilds records with raw features only."""
    import types
    r = types.SimpleNamespace(dicer_call="24", frac_top=0.5, pct_1u=0.1,
                              pct_5a=0.8, pingpong_frac=0.0, phase_score=0.0,
                              reads=500)
    RuleScorer(CFG).decide(r)
    assert r.score == 0.0
