"""Both strands producing is the siRNA signature -- not a 50/50 split.

A cluster can be genuinely bidirectional at 70/30 (transcription and AGO loading
have their own biases). The equivalence test asked for symmetry instead, and was
underpowered on top of it: a perfect 50/50 locus with 50 distinct sequences
scored 1.68, below the threshold of 2 it was meant to pass.
"""
from evidence import strand_bidirectional_score


def test_balanced_locus_scores_high():
    assert strand_bidirectional_score(0.5, 100) > 5


def test_asymmetric_but_bidirectional_still_passes():
    # 70/30 with plenty of reads: both strands clearly produce
    assert strand_bidirectional_score(0.7, 100) > 2


def test_direction_does_not_matter():
    assert (strand_bidirectional_score(0.7, 100)
            == strand_bidirectional_score(0.3, 100))


def test_uni_strand_locus_scores_zero():
    assert strand_bidirectional_score(1.0, 100) == 0.0


def test_nearly_uni_strand_stays_low():
    # 2% on the minority strand is below the 10% null: no evidence
    assert strand_bidirectional_score(0.98, 100) < 2


def test_it_has_power_where_the_equivalence_test_did_not():
    from evidence import strand_balance_score
    assert strand_balance_score(0.5, 50) < 2          # the old term fails here
    assert strand_bidirectional_score(0.5, 50) > 2    # this one does not


def test_small_n_gives_no_evidence():
    assert strand_bidirectional_score(0.5, 4) < 2


def test_guards():
    assert strand_bidirectional_score(None, 100) == 0.0
    assert strand_bidirectional_score(0.5, 0) == 0.0
