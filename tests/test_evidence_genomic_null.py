"""The 5' base null is the locus's own DNA composition, not 0.25.

A cut placed at random in the locus yields a 5'U as often as the locus has T on
the read's strand. Assuming 0.25 ignores that: a T-rich transposon would look
1U-biased with no selection at all. Measured, the genomic null is 0.272 in
Drosophila and 0.332 in Arabidopsis -- and reads' own positions 3-5, which look
like a cheaper proxy, disagree with the DNA in plants (0.255 against 0.332).

The null is strand-weighted because pct_1u pools reads from both strands: a read
on the minus strand starts with T where the plus strand carries A.
"""
import pytest

from samio import base_null


def test_all_T_on_the_plus_strand():
    # every plus-strand read starts with T; every minus-strand read with A
    assert base_null("TTTTTTTTTT", "T", frac_top=1.0) == 1.0
    assert base_null("TTTTTTTTTT", "T", frac_top=0.0) == 0.0


def test_balanced_strands_average_the_two():
    assert base_null("TTTTTTTTTT", "T", frac_top=0.5) == 0.5


def test_composition_is_counted_exactly():
    # ACGTACGTAC holds 2 T in 10 bases
    assert base_null("ACGTACGTAC", "T", frac_top=1.0) == pytest.approx(0.2, abs=0.001)


def test_a_T_rich_locus_raises_the_bar():
    """This is the point: 40% T means 40% of random cuts start with U."""
    seq = "T" * 40 + "A" * 20 + "C" * 20 + "G" * 20
    assert base_null(seq, "T", frac_top=1.0) == pytest.approx(0.40, abs=0.01)


def test_A_null_mirrors_T_on_the_other_strand():
    seq = "T" * 40 + "A" * 20 + "C" * 20 + "G" * 20
    assert base_null(seq, "A", frac_top=0.0) == pytest.approx(0.40, abs=0.01)


def test_lowercase_counts_and_N_is_excluded():
    """N cannot be cut into a 5' base, so it leaves the denominator."""
    assert base_null("ttttNNNN", "T", frac_top=1.0) == pytest.approx(1.0, abs=0.001)
    assert base_null("ttttcccc", "T", frac_top=1.0) == pytest.approx(0.5, abs=0.001)


def test_empty_sequence_falls_back_to_uniform():
    assert base_null("", "T", frac_top=1.0) == 0.25
