"""The plant config must separate the three canonical plant classes.

Not a validation of the weights (no plant dataset has been run yet) — these
tests pin the intended behaviour so a future calibration cannot silently break
the discrimination the config exists for.
"""
from records import ClusterRecord
from decider import RuleScorer, load_config


def _cfg():
    return load_config("rules.plant.yaml")


def _decide(**kw):
    base = dict(chrom="c", start=0, end=1000, reads=500, frac_top=0.5,
                pct_1u=0.2, pct_5a=0.2, pingpong_frac=0.0, phase_score=0.0)
    base.update(kw)
    r = ClusterRecord(**base)
    RuleScorer(_cfg()).decide(r)
    return r


def test_phased_21nt_with_5u_is_phasirna():
    r = _decide(dicer_call="21", phase_score=20.0, pct_1u=0.9, pct_5a=0.05)
    assert r.cls == "phasiRNA"


def test_unphased_21nt_dualstrand_is_plain_sirna():
    r = _decide(dicer_call="21", phase_score=0.0, pct_1u=0.3, frac_top=0.5)
    assert r.cls == "siRNA"


def test_24nt_with_5a_is_hcsirna():
    r = _decide(dicer_call="24", pct_5a=0.9, pct_1u=0.05)
    assert r.cls == "hcsiRNA"


def test_no_pirna_class_exists():
    # plants have no PIWI: a 27-nt locus must not be called piRNA, because the
    # class is not in the config at all
    assert "piRNA" not in _cfg()["classes"]
    r = _decide(dicer_call="27", pct_1u=0.9)
    assert r.cls in ("phasiRNA", "hcsiRNA", "siRNA")


def test_ping_pong_is_not_scored_in_plants():
    # no PIWI, no ping-pong cycle: the signal must not move any plant score
    quiet = _decide(dicer_call="24", pct_5a=0.9, pingpong_frac=0.0)
    loud = _decide(dicer_call="24", pct_5a=0.9, pingpong_frac=0.9)
    assert quiet.class_scores == loud.class_scores
