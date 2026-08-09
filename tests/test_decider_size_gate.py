"""Every evidence term is conditioned on the size fitting the class.

Ping-pong and phasing were gated by the class's own size window; 5'U, 5'A and
the dual-strand term were not. Measured consequence, per class, of loci called
into a class whose size window they sit outside:

    Drosophila  piRNA 14.5% (carried by w_1u)      siRNA 16.3% (w_fractop_mid)
    Tribolium   siRNA 22.0%                        piRNA  1.6%
    Arabidopsis phasiRNA 20.7%                     hcsiRNA 2.9% (w_5a)
    rice        siRNA 27.3%   hcsiRNA 12.1%        phasiRNA 13.0%
    rice        phasiRNA24 0.5%  <- the one class whose terms were all gated

In the worst group, 675 Drosophila loci of 18-20 nt, the size and ping-pong
terms are both 0 and the piRNA call rests entirely on a 5'U of 0.78 -- higher
than the 0.70 of real piRNA, because a 5'U is common to many small RNAs and
only means piRNA together with the size. A 19-nt RNA is not a piRNA: the window
comes from what Piwi loads.
"""
import pytest

from decider import RuleScorer
from records import ClusterRecord


CFG = {
    "classes": {
        "piRNA": {"size": {"lo": 24, "hi": 31}, "w_size": 1.0, "w_1u": 1.0,
                  "w_pingpong": 1.5, "defining_evidence": [["1u", "pingpong"]]},
        "siRNA": {"size": {"lo": 20, "hi": 22}, "w_size": 1.0,
                  "w_fractop_mid": 1.0, "defining_evidence": ["fractop_bid"]},
    },
    "decision": {"review_confidence": 0.5, "min_margin": 0.1, "min_reads": 10},
}


def _rec(dicer, pct_1u=0.0, frac_top=0.5, pct_5a=0.0):
    return ClusterRecord(chrom="c", start=0, end=400, reads=500,
                         dicer_call=str(dicer), pct_1u=pct_1u,
                         frac_top=frac_top, pct_5a=pct_5a, pingpong_frac=0.0)


def _term(scorer, cls, record, name):
    spec = scorer.cfg["classes"][cls]
    return scorer._score_class(cls, spec, record)[1].get(name, 0.0)


def test_the_5u_term_is_zero_where_the_size_does_not_fit():
    """19 nt is 5 below the piRNA window: no size support, no 5'U credit."""
    scorer = RuleScorer(CFG)
    assert _term(scorer, "piRNA", _rec(19, pct_1u=0.9), "1u") == 0.0


def test_the_5u_term_survives_inside_the_window():
    scorer = RuleScorer(CFG)
    assert _term(scorer, "piRNA", _rec(27, pct_1u=0.9), "1u") == pytest.approx(0.9)


def test_the_5u_term_decays_at_the_edge():
    """23 nt is one below the window: partial credit, not full and not none."""
    scorer = RuleScorer(CFG)
    edge = _term(scorer, "piRNA", _rec(23, pct_1u=0.9), "1u")
    assert 0.0 < edge < 0.9


def test_the_dual_strand_term_is_gated_too():
    scorer = RuleScorer(CFG)
    inside = _term(scorer, "siRNA", _rec(21, frac_top=0.5), "fractop_mid")
    outside = _term(scorer, "siRNA", _rec(30, frac_top=0.5), "fractop_mid")
    assert inside > 0.0
    assert outside == 0.0


def test_the_5a_term_is_gated_too():
    cfg = {"classes": {"hcsiRNA": {"size": {"lo": 23, "hi": 24}, "w_size": 1.0,
                                   "w_5a": 1.0, "defining_evidence": ["5a"]}},
           "decision": {"review_confidence": 0.5, "min_margin": 0.1,
                        "min_reads": 10}}
    scorer = RuleScorer(cfg)
    assert _term(scorer, "hcsiRNA", _rec(24, pct_5a=0.8), "5a") == pytest.approx(0.8)
    assert _term(scorer, "hcsiRNA", _rec(18, pct_5a=0.8), "5a") == 0.0


def test_a_19nt_locus_with_a_high_5u_is_not_a_pirna():
    """The 675-locus population that motivated this, as one record."""
    scorer = RuleScorer(CFG)
    r = _rec(19, pct_1u=0.78, frac_top=0.68)
    scorer.decide(r)
    assert r.cls != "piRNA"


def test_a_real_pirna_is_untouched():
    scorer = RuleScorer(CFG)
    r = _rec(27, pct_1u=0.70, frac_top=0.56)
    scorer.decide(r)
    assert r.cls == "piRNA"


def test_a_locus_with_no_size_call_scores_nothing_anywhere():
    """dicer_call 'N' already zeroes the score; the class terms follow."""
    scorer = RuleScorer(CFG)
    r = _rec("N", pct_1u=0.9)
    assert _term(scorer, "piRNA", r, "1u") == 0.0
