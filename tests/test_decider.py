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

def test_a_size_between_the_windows_is_ambiguous_however_strong_the_evidence():
    """23 nt sits one below piRNA and one above siRNA, so neither holds.

    This used to assert that ping-pong picked the class here. It no longer can:
    a class is asserted only over a size inside its own window, and no amount of
    evidence changes which sizes a class is defined over.
    """
    scorer = RuleScorer(_cfg())
    base = dict(chrom="c", start=0, end=100, reads=500, dicer_call="23",
                frac_top=0.5, pct_1u=0.5)
    no_pp = ClusterRecord(**base, pingpong_frac=0.0)
    strong_pp = ClusterRecord(**base, pingpong_frac=0.5)
    scorer.decide(no_pp)
    scorer.decide(strong_pp)
    assert no_pp.cls == "ambiguous"
    assert strong_pp.cls == "ambiguous"
    # the trail is kept: ping-pong still moved which class came closest
    assert "nearest=siRNA" in no_pp.evidence
    assert "nearest=piRNA" in strong_pp.evidence

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


def test_no_size_call_scores_zero():
    # dicer_call == "N" means no size call, so no class evidence is meaningful,
    # even if a class is still picked and non-zero defining evidence exists.
    r = ClusterRecord(chrom="c", start=0, end=100, reads=500,
                      dicer_call="N", frac_top=0.95, pct_1u=0.9, pingpong_frac=0.3)
    r.e_1u, r.e_pingpong = 0.9, 0.9
    RuleScorer(_cfg()).decide(r)
    # No size means no class window contains it: reported ambiguous, and a
    # locus with nothing asserted over it carries no score rather than a zero.
    assert r.cls == "ambiguous"
    assert r.score is None


def test_pirna_defining_evidence_is_1u_or_pingpong():
    # piRNA is size + (5'U OR ping-pong): either evidence alone must carry the score.
    r_1u = ClusterRecord(chrom="c", start=0, end=100, reads=500, dicer_call="27")
    r_1u.e_1u, r_1u.e_pingpong = 0.8, 0.0
    r_pp = ClusterRecord(chrom="c", start=0, end=100, reads=500, dicer_call="27")
    r_pp.e_1u, r_pp.e_pingpong = 0.0, 0.8
    cfg = _cfg()
    spec = cfg["classes"]["piRNA"]
    # An instance now: `size` evidence is scored against the class, so the
    # lookup needs the config, not just the record.
    scorer = RuleScorer(cfg)
    # The OR slot is the first entry; `size` follows it and scores 0 here,
    # since these bare records carry no size_counts.
    assert scorer.defining_scores(spec, r_1u)[0] == 0.8
    assert scorer.defining_scores(spec, r_pp)[0] == 0.8
