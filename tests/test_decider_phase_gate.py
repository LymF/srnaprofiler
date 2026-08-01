"""Phasing is size-gated, which is what lets one config serve plants and arthropods.

Ungated, w_phase was identical in both classes, so it added the same constant to
both scores and changed neither the argmax nor the two-way softmax: phasing was
inert. Gated, a phased 21-nt locus reinforces siRNA (plant phasiRNA/tasiRNA) and a
phased 27-nt locus reinforces piRNA (arthropod Zucchini products).
"""
from records import ClusterRecord
from decider import RuleScorer, load_config


def _cfg():
    return load_config("rules.yaml")


def _score(dicer_call, phase_score):
    r = ClusterRecord(chrom="c", start=0, end=1000, reads=500, dicer_call=dicer_call,
                      frac_top=0.5, pct_1u=0.5, pingpong_frac=0.0,
                      phase_score=phase_score)
    RuleScorer(_cfg()).decide(r)
    return r


def test_phasing_reinforces_sirna_at_21nt():
    flat = _score("21", 0.0)
    phased = _score("21", 20.0)
    gain_si = phased.class_scores["siRNA"] - flat.class_scores["siRNA"]
    gain_pi = phased.class_scores["piRNA"] - flat.class_scores["piRNA"]
    assert gain_si > 0.4          # inside the siRNA window: full weight
    assert gain_pi == 0.0         # far outside the piRNA window: gated off
    assert phased.cls == "siRNA"


def test_phasing_reinforces_pirna_at_27nt():
    flat = _score("27", 0.0)
    phased = _score("27", 20.0)
    gain_pi = phased.class_scores["piRNA"] - flat.class_scores["piRNA"]
    gain_si = phased.class_scores["siRNA"] - flat.class_scores["siRNA"]
    assert gain_pi > 0.4
    assert gain_si == 0.0
    assert phased.cls == "piRNA"
