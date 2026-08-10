"""A class is only asserted when the size falls inside its window.

The scorer used to pick the nearest class however far the size sat from it: a
19-nt locus became piRNA, a 23-nt one too. `review` could not carry that, because
review means "the classes compete" and this is "none of them fits" -- two
different statements the roadmap has kept apart since 2026-08-07.

The cost is known and named: two of the 90 reference piRNA clusters are
represented by a single locus of 78 and 46 reads whose sizes spread from 18 to
27 nt with a 19% and 23% mode at 23. Under this rule they are reported
ambiguous, so piC recall reads 83/90 instead of 85/90. The loci do not vanish --
they are reported, without a class asserted over them. 85/90 with two forced
calls is a worse claim than 83/90 with two honest ones.

Purity was measured as an alternative and rejected: sRNAminer's minRepeatRatio
of 0.5, applied here, drops piC recall to 66/90, because size purity is a
property of the kingdom -- median 0.853 for Arabidopsis hcsiRNA against 0.410
for Drosophila piRNA. A piRNA locus spreads over 24-31 nt by construction;
demanding purity of it is demanding that it be an siRNA.
"""
import pytest

from decider import RuleScorer
from records import ClusterRecord


CFG = {
    "classes": {
        "piRNA": {"size": {"lo": 24, "hi": 31}, "w_size": 1.0, "w_1u": 1.0,
                  "defining_evidence": [["1u", "pingpong"]]},
        "siRNA": {"size": {"lo": 20, "hi": 22}, "w_size": 1.0,
                  "w_fractop_mid": 1.0, "defining_evidence": ["fractop_bid"]},
    },
    "decision": {"review_confidence": 0.5, "min_margin": 0.1, "min_reads": 10},
}


def _rec(dicer, pct_1u=0.8, frac_top=0.5, reads=500):
    return ClusterRecord(chrom="c", start=0, end=400, reads=reads,
                         dicer_call=str(dicer), pct_1u=pct_1u,
                         frac_top=frac_top, pingpong_frac=0.0)


def test_a_size_inside_the_window_keeps_its_class():
    r = _rec(27)
    RuleScorer(CFG).decide(r)
    assert r.cls == "piRNA"


def test_one_nt_outside_is_already_ambiguous():
    """23 nt sits one below piRNA and one above siRNA: neither definition holds."""
    r = _rec(23)
    RuleScorer(CFG).decide(r)
    assert r.cls == "ambiguous"


def test_far_outside_is_ambiguous():
    r = _rec(19)
    RuleScorer(CFG).decide(r)
    assert r.cls == "ambiguous"


def test_an_ambiguous_locus_reports_no_score():
    """Nothing was asserted, so there is no evidence strength to report."""
    r = _rec(19)
    RuleScorer(CFG).decide(r)
    assert r.score is None


def test_an_ambiguous_locus_is_flagged_for_review():
    """review means "do not trust this call". An ambiguous one qualifies.

    Clearing it would be worse than redundant: anyone filtering review == False
    for the confident calls would collect the ambiguous loci along with them.
    """
    r = _rec(19)
    RuleScorer(CFG).decide(r)
    assert r.review is True


def test_the_nearest_class_is_kept_in_the_evidence():
    """The trail must survive: which class came closest, and why."""
    r = _rec(19)
    RuleScorer(CFG).decide(r)
    assert "nearest=" in r.evidence


def test_a_locus_with_no_size_call_is_ambiguous():
    """With every term size-gated, an 'N' scores zero for every class, so the
    argmax would be whichever class the config happens to list first."""
    r = _rec("N")
    RuleScorer(CFG).decide(r)
    assert r.cls == "ambiguous"
    assert r.score is None


def test_the_class_scores_are_still_reported():
    """Ambiguous is a verdict about the size, not a refusal to measure."""
    r = _rec(19)
    RuleScorer(CFG).decide(r)
    assert set(r.class_scores) == {"piRNA", "siRNA"}


def test_a_class_without_a_size_window_is_never_ambiguous():
    cfg = {"classes": {"anything": {"w_1u": 1.0,
                                    "defining_evidence": ["1u"]}},
           "decision": {"review_confidence": 0.0, "min_margin": 0.0,
                        "min_reads": 10}}
    r = _rec(19)
    RuleScorer(cfg).decide(r)
    assert r.cls == "anything"
