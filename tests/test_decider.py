from records import ClusterRecord
from decider import RuleScorer, load_config

def _cfg():
    return load_config("rules.yaml")

def test_pirna_strong_signal():
    r = ClusterRecord(chrom="c", start=0, end=100, reads=500,
                      dicer_call="27", frac_top=0.95, pct_1u=0.9, pingpong_z=8.0)
    RuleScorer(_cfg()).decide(r)
    assert r.cls == "piRNA"
    assert r.confidence > 0.5
    assert r.review is False

def test_sirna_dualstrand_21nt():
    r = ClusterRecord(chrom="c", start=0, end=100, reads=500,
                      dicer_call="21", frac_top=0.5, pct_1u=0.2, pingpong_z=0.0)
    RuleScorer(_cfg()).decide(r)
    assert r.cls == "siRNA"

def test_weak_locus_flagged_review():
    r = ClusterRecord(chrom="c", start=0, end=100, reads=5,
                      dicer_call="N", frac_top=0.5)
    RuleScorer(_cfg()).decide(r)
    assert r.review is True
