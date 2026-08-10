"""A class can be switched off in the config, and then it is off everywhere.

phasiRNA24 describes a monocot pathway: 24-nt reproductive phasiRNA come from
DCL5, which rice and maize have and Arabidopsis does not. Shipping it enabled in
the generic plant config gave 3514 loci and 46.4% of the read mass in a species
that lacks the pathway -- loci with a median phase_pvalue of 1.000, which are
simply the deepest hc-siRNA loci.

Switching a class off has to reach the clustering too. Class size windows feed
`clustering._same_size_group`, which decides whether two islands may merge, so
merely skipping the class in the scorer would leave the loci themselves changed:
adding the {24,24} window moved the Arabidopsis run from 28,754 to 32,992 loci.
"""
import pytest

from classify import class_size_windows
from decider import RuleScorer
from records import ClusterRecord


def _cfg(enabled=None):
    spec = {"size": {"lo": 24, "hi": 24}, "w_size": 1.0, "w_phase": 4.0,
            "phase_on": "cycles", "defining_evidence": ["phase", "size"]}
    if enabled is not None:
        spec["enabled"] = enabled
    return {
        "classes": {
            "hcsiRNA": {"size": {"lo": 23, "hi": 24}, "w_size": 1.0,
                        "w_5a": 1.0, "defining_evidence": ["5a"]},
            "phasiRNA24": spec,
        },
        "decision": {"review_confidence": 0.5, "min_margin": 0.1,
                     "min_reads": 10},
    }


def _rec(number=11, pct_5a=0.5):
    r = ClusterRecord(chrom="c", start=0, end=400, reads=500, dicer_call="24",
                      pct_5a=pct_5a, frac_top=0.5)
    r.phase_number = number
    return r


def test_a_disabled_class_is_never_assigned():
    r = _rec()
    RuleScorer(_cfg(enabled=False)).decide(r)
    assert r.cls == "hcsiRNA"


def test_the_same_locus_takes_the_class_when_enabled():
    """Without this the first test would pass for the wrong reason."""
    r = _rec()
    RuleScorer(_cfg(enabled=True)).decide(r)
    assert r.cls == "phasiRNA24"


def test_a_class_with_no_flag_is_enabled():
    """Existing configs carry no `enabled` key and must not change."""
    r = _rec()
    RuleScorer(_cfg()).decide(r)
    assert r.cls == "phasiRNA24"


def test_a_disabled_class_leaves_no_score_behind():
    r = _rec()
    RuleScorer(_cfg(enabled=False)).decide(r)
    assert "phasiRNA24" not in r.class_scores


def test_a_disabled_class_does_not_shape_the_clustering():
    """Its size window must not reach _same_size_group, or the loci differ."""
    assert (24, 24) not in class_size_windows(_cfg(enabled=False))
    assert (24, 24) in class_size_windows(_cfg(enabled=True))


def test_disabling_every_class_is_refused():
    """An empty class set would make every locus ambiguous, silently."""
    cfg = _cfg(enabled=False)
    cfg["classes"]["hcsiRNA"]["enabled"] = False
    with pytest.raises(ValueError, match="no enabled class"):
        RuleScorer(cfg)
