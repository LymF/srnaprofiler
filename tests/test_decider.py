import pytest

from records import ClusterRecord
from decider import RuleScorer, load_config

def _cfg():
    return load_config("rules.yaml")

def test_pirna_strong_signal():
    r = ClusterRecord(chrom="c", start=0, end=100, reads=500,
                      dicer_call="27", frac_top=0.95, pct_1u=0.9, pingpong_frac=0.3)
    RuleScorer(_cfg()).decide(r)
    assert r.cls == "piRNA"
    assert r.confidence > 0.5
    assert r.review is False

def test_sirna_dualstrand_21nt():
    r = ClusterRecord(chrom="c", start=0, end=100, reads=500,
                      dicer_call="21", frac_top=0.5, pct_1u=0.2, pingpong_frac=0.0)
    RuleScorer(_cfg()).decide(r)
    assert r.cls == "siRNA"

def test_pingpong_frac_flips_ambiguous_to_pirna():
    # a 23-nt dual-strand locus sits between the size windows; ping-pong decides.
    scorer = RuleScorer(_cfg())
    base = dict(chrom="c", start=0, end=100, reads=500, dicer_call="23",
                frac_top=0.5, pct_1u=0.5)
    no_pp = ClusterRecord(**base, pingpong_frac=0.0)
    strong_pp = ClusterRecord(**base, pingpong_frac=0.5)
    scorer.decide(no_pp)
    scorer.decide(strong_pp)
    assert no_pp.cls == "siRNA"       # without ping-pong -> siRNA
    assert strong_pp.cls == "piRNA"   # strong ping-pong -> piRNA

def test_pingpong_does_not_override_sirna_size():
    # 21-nt dual-strand loci (canonical endo-siRNA) pick up 10-nt 5'-5' overlaps by
    # geometry alone, so ping-pong must not drag them into piRNA: the size gate zeroes
    # both the piRNA bonus and the siRNA penalty outside each class's size window.
    r = ClusterRecord(chrom="c", start=0, end=100, reads=500, dicer_call="21",
                      frac_top=0.5, pct_1u=0.65, pingpong_frac=0.5)
    RuleScorer(_cfg()).decide(r)
    assert r.cls == "siRNA"


def test_size_gate_leaves_real_pirna_untouched():
    # A 27-nt locus is inside the piRNA window, so the gate is 1.0 and ping-pong
    # keeps its full weight.
    scorer = RuleScorer(_cfg())
    base = dict(chrom="c", start=0, end=100, reads=500, dicer_call="27",
                frac_top=0.5, pct_1u=0.5)
    no_pp = ClusterRecord(**base, pingpong_frac=0.0)
    strong_pp = ClusterRecord(**base, pingpong_frac=0.5)
    scorer.decide(no_pp)
    scorer.decide(strong_pp)
    assert strong_pp.cls == "piRNA"
    assert strong_pp.class_scores["piRNA"] - no_pp.class_scores["piRNA"] == pytest.approx(1.5)


def test_weak_locus_flagged_review():
    r = ClusterRecord(chrom="c", start=0, end=100, reads=5,
                      dicer_call="N", frac_top=0.5)
    RuleScorer(_cfg()).decide(r)
    assert r.review is True
