"""Size is what decides the class, but it never entered the score.

`w_size` multiplies the class window in `_score_class`, so the size call is what
picks the class -- but `defining_evidence` never lists it, so `score`, the column
that says how sure the call is, ignores the size entirely. A 24-nt locus with a
perfect size profile and a weak 5'A reports the same 0.0 as a locus with no
evidence at all.

e_size is the first class-dependent strength: the other e_* are computed once in
features, but a size only means something against a given class's window, so it
is computed in the scorer.

  p_hat = mass inside the class window / mass inside the extractor window
  n     = distinct_seqs, the same denominator the other strengths use
  p0    = width of the class window / width of the extractor window,
          i.e. what a locus with no size preference would land there by chance
"""
import pytest

from decider import RuleScorer
from records import ClusterRecord


CFG = {
    "classes": {
        "hcsiRNA": {"size": {"lo": 23, "hi": 24}, "w_size": 1.0,
                    "defining_evidence": ["size"]},
    },
    "decision": {"review_confidence": 0.5, "min_margin": 0.1, "min_reads": 10},
    "extractor": {"dicermin": 18, "dicermax": 35},
}


def _rec(counts, distinct=100, dicer="24"):
    r = ClusterRecord(chrom="c", start=0, end=500)
    r.size_counts = dict(counts)
    r.distinct_seqs = distinct
    r.dicer_call = dicer
    return r


def test_a_locus_entirely_in_the_window_has_full_size_evidence():
    scorer = RuleScorer(CFG)
    r = _rec({24: 100.0})
    assert scorer.size_strength(CFG["classes"]["hcsiRNA"], r) > 0.9


def test_a_locus_outside_the_window_has_none():
    scorer = RuleScorer(CFG)
    r = _rec({21: 100.0}, dicer="21")
    assert scorer.size_strength(CFG["classes"]["hcsiRNA"], r) == 0.0


def test_a_mixed_locus_lands_in_between():
    scorer = RuleScorer(CFG)
    pure = scorer.size_strength(CFG["classes"]["hcsiRNA"], _rec({24: 100.0}))
    mixed = scorer.size_strength(CFG["classes"]["hcsiRNA"],
                                 _rec({24: 50.0, 21: 50.0}))
    assert 0.0 < mixed < pure


def test_the_null_is_the_share_of_the_window_a_random_size_would_hit():
    """23-24 is 2 of the 18 lengths the extractor keeps, so p0 = 2/18."""
    scorer = RuleScorer(CFG)
    # A locus at exactly the null strength must score 0.
    counts = {L: 100.0 for L in range(18, 36)}
    assert scorer.size_strength(CFG["classes"]["hcsiRNA"], _rec(counts)) == 0.0


def test_a_small_locus_is_penalised_for_its_own_uncertainty():
    """Wilson lower bound: the same fraction with fewer sequences scores less."""
    scorer = RuleScorer(CFG)
    many = scorer.size_strength(CFG["classes"]["hcsiRNA"],
                                _rec({24: 100.0}, distinct=200))
    few = scorer.size_strength(CFG["classes"]["hcsiRNA"],
                               _rec({24: 100.0}, distinct=5))
    assert few < many


def test_an_empty_locus_scores_zero():
    scorer = RuleScorer(CFG)
    assert scorer.size_strength(CFG["classes"]["hcsiRNA"], _rec({}, distinct=0)) == 0.0


def test_size_is_accepted_in_defining_evidence():
    """The guard refuses unknown names; 'size' must not be one of them."""
    RuleScorer(CFG)


def test_the_decision_reads_size_from_defining_evidence():
    """End to end: a locus of the right size gets a score, not a zero."""
    scorer = RuleScorer(CFG)
    r = _rec({24: 100.0})
    r.reads = 100
    scorer.decide(r)
    assert r.cls == "hcsiRNA"
    assert r.score > 0.9


def test_a_class_without_size_in_its_evidence_is_unaffected():
    cfg = {
        "classes": {"siRNA": {"size": {"lo": 21, "hi": 22}, "w_size": 1.0,
                              "defining_evidence": ["fractop_bid"]}},
        "decision": {"review_confidence": 0.5, "min_margin": 0.1, "min_reads": 10},
    }
    scorer = RuleScorer(cfg)
    r = _rec({21: 100.0}, dicer="21")
    r.reads = 100
    r.e_fractop_bid = 0.42
    scorer.decide(r)
    assert r.score == 0.42
