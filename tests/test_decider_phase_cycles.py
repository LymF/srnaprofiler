"""The class term for phasing can read cycles instead of abundance.

`w_phase * sigmoid(phase_score)` is a depth detector, measured: adding a 24-nt
phasiRNA class to the plant config moved 2006 hc-siRNA loci in Arabidopsis --
which has no DCL5 and therefore no 24-nt phasiRNA at all. Those loci are 5x
deeper than the ones that stayed (555 vs 112 reads), their phase_pvalue is
1.000, and they hit p <= 0.01 at 1.1% against a 0.9% background. They are deep,
not phased.

The count of cycles in register is the quantity that discriminates (AUC 0.826
against size-matched controls, versus 0.560 for the p-value), so a class can ask
for that instead with `phase_on: cycles`. The default stays `score`, so no
existing class changes.
"""
import pytest

from decider import RuleScorer
from records import ClusterRecord


def _cfg(phase_on=None, w_phase=3.0):
    spec = {"size": {"lo": 24, "hi": 24}, "w_size": 1.0, "w_phase": w_phase,
            "defining_evidence": ["phase"]}
    if phase_on:
        spec["phase_on"] = phase_on
    return {"classes": {"phasiRNA24": spec, "hcsiRNA": {
                "size": {"lo": 23, "hi": 24}, "w_size": 1.0, "w_5a": 1.0,
                "defining_evidence": ["5a"]}},
            "decision": {"review_confidence": 0.5, "min_margin": 0.1,
                         "min_reads": 10}}


def _rec(number=0, score=0.0, dicer="24", reads=500):
    r = ClusterRecord(chrom="c", start=0, end=500, reads=reads,
                      dicer_call=dicer, pct_5a=0.5)
    r.phase_number, r.phase_score = number, score
    return r


def test_cycles_mode_ignores_the_abundance_driven_score():
    """Same cycles, wildly different phase_score -> same phase term."""
    scorer = RuleScorer(_cfg(phase_on="cycles"))
    spec = scorer.cfg["classes"]["phasiRNA24"]
    shallow = scorer._score_class("phasiRNA24", spec, _rec(number=8, score=2.0))
    deep = scorer._score_class("phasiRNA24", spec, _rec(number=8, score=20.0))
    assert shallow[0] == deep[0]


def test_cycles_mode_grows_with_the_run_in_register():
    scorer = RuleScorer(_cfg(phase_on="cycles"))
    spec = scorer.cfg["classes"]["phasiRNA24"]
    few = scorer._score_class("phasiRNA24", spec, _rec(number=2))[0]
    many = scorer._score_class("phasiRNA24", spec, _rec(number=10))[0]
    assert many > few


def test_a_deep_unphased_locus_gets_no_phase_credit_in_cycles_mode():
    """The exact profile of the 2006: high phase_score, two cycles."""
    scorer = RuleScorer(_cfg(phase_on="cycles"))
    spec = scorer.cfg["classes"]["phasiRNA24"]
    _, terms = scorer._score_class("phasiRNA24", spec,
                                   _rec(number=2, score=8.4, reads=555))
    assert terms["phase"] < 0.1


def test_the_default_still_reads_phase_score():
    """No `phase_on` means the shipped behaviour, unchanged."""
    scorer = RuleScorer(_cfg())
    spec = scorer.cfg["classes"]["phasiRNA24"]
    low = scorer._score_class("phasiRNA24", spec, _rec(number=2, score=1.0))[0]
    high = scorer._score_class("phasiRNA24", spec, _rec(number=2, score=20.0))[0]
    assert high > low


def test_an_unknown_phase_mode_is_refused_at_construction():
    """A typo must not silently fall back to the depth detector."""
    with pytest.raises(ValueError, match="phase_on"):
        RuleScorer(_cfg(phase_on="abundance"))


def test_cycles_mode_is_still_size_gated():
    """A 21-nt locus must not collect the 24-nt class's phasing credit."""
    scorer = RuleScorer(_cfg(phase_on="cycles"))
    spec = scorer.cfg["classes"]["phasiRNA24"]
    right = scorer._score_class("phasiRNA24", spec, _rec(number=10, dicer="24"))[0]
    wrong = scorer._score_class("phasiRNA24", spec, _rec(number=10, dicer="21"))[0]
    assert wrong < right


def test_a_locus_with_no_size_call_gets_no_cycle_credit():
    scorer = RuleScorer(_cfg(phase_on="cycles"))
    spec = scorer.cfg["classes"]["phasiRNA24"]
    _, terms = scorer._score_class("phasiRNA24", spec, _rec(number=10, dicer="N"))
    assert terms.get("phase", 0.0) == 0.0


def test_a_genuinely_phased_24nt_locus_beats_hcsirna():
    """Specificity is worthless if the class can never win.

    Arabidopsis has no DCL5, so zero 24-nt phasiRNA there is the right answer --
    but a locus with a full run in register must still take the class, or the
    rule is dead weight rather than a specific one.
    """
    scorer = RuleScorer(_cfg(phase_on="cycles", w_phase=4.0))
    phased = _rec(number=11, score=3.0)
    scorer.decide(phased)
    assert phased.cls == "phasiRNA24"

    unphased = _rec(number=2, score=12.0)   # deep, not phased
    scorer.decide(unphased)
    assert unphased.cls == "hcsiRNA"
