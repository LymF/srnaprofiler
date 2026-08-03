"""Null tests behind the per-locus evidence scores.

score = log10(1/p), so score 2 means p = 0.01. p0 = 0.25 for the 5' base
follows proTRAC, which states it tests 1T/10A "based on random base
composition (25% for each nucleotide)".
"""
import math

from evidence import score_from_p, base_bias_score, strand_bias_score


def test_score_from_p_is_log10_of_one_over_p():
    assert score_from_p(0.01) == 2.0
    assert score_from_p(0.001) == 3.0
    assert score_from_p(1.0) == 0.0


def test_score_from_p_never_negative_and_caps_at_zero_p():
    assert score_from_p(1.5) == 0.0          # p > 1 from rounding
    assert score_from_p(0.0) == 300.0        # underflow, capped


def test_base_bias_strong_signal_scores_high():
    # 90 of 100 reads starting with the same base, against 25% by chance
    assert base_bias_score(90, 100) > 20


def test_base_bias_at_chance_scores_zero():
    assert base_bias_score(25, 100) < 0.5


def test_base_bias_shallow_locus_cannot_prove_much():
    # same 90% proportion, but only 10 reads
    deep = base_bias_score(90, 100)
    shallow = base_bias_score(9, 10)
    assert shallow < deep
    assert shallow < 6


def test_base_bias_handles_degenerate_input():
    assert base_bias_score(0, 0) == 0.0
    assert base_bias_score(5, 3) == 0.0      # k > n is nonsense


def test_strand_bias_scores_asymmetry_either_way():
    assert strand_bias_score(0.95, 200) > 10
    assert strand_bias_score(0.05, 200) > 10


def test_strand_bias_zero_for_balanced_locus():
    assert strand_bias_score(0.5, 200) < 0.5
