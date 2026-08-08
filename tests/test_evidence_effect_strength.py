"""Strength of an effect, not its significance.

A p-value mixes effect size with sample size: with a big n, a tiny bias gives
an astronomical p. Measured on the piRNA loci, `score` correlated 0.629 with
distinct_seqs and only 0.323 with the effect, and loci scoring 50+ sat at 24.9%
inside known piRNA clusters, below the 28.3% background.

This returns the Wilson lower bound of k/n, rescaled by how much room the null
leaves above it, so 0.61 reads as "at least 61% of the way from the null to
certainty, at 95% confidence".
"""
import pytest

from evidence import effect_strength


def test_the_same_effect_scores_alike_at_different_n():
    """The property the p-value score does not have."""
    small = effect_strength(38, 50, 0.25)
    large = effect_strength(760, 1000, 0.25)
    assert abs(large - small) < 0.25


def test_more_data_never_lowers_the_strength():
    """A tighter interval can only move the lower bound up."""
    assert effect_strength(760, 1000, 0.25) >= effect_strength(38, 50, 0.25)


def test_effect_at_the_null_is_zero():
    assert effect_strength(250, 1000, 0.25) == 0.0


def test_effect_below_the_null_is_zero():
    assert effect_strength(100, 1000, 0.25) == 0.0


def test_a_perfect_effect_with_plenty_of_data_approaches_one():
    assert effect_strength(2000, 2000, 0.25) > 0.95


def test_a_perfect_effect_with_little_data_does_not():
    """Three out of three is not evidence of much."""
    assert effect_strength(3, 3, 0.25) < 0.7


def test_the_null_rescales_the_range():
    """Same observed proportion is worth more against a harder null."""
    assert effect_strength(600, 1000, 0.10) > effect_strength(600, 1000, 0.50)


def test_guards():
    assert effect_strength(10, 0, 0.25) == 0.0
    assert effect_strength(-1, 100, 0.25) == 0.0
    assert effect_strength(101, 100, 0.25) == 0.0
    assert effect_strength(50, 100, 1.0) == 0.0
    assert effect_strength(50, 100, None) == 0.0
