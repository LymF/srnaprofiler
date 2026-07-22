import math
from signals import phase_metrics, hypergeometric_phase_pvalue

def test_phase_metrics_perfectly_phased():
    # three reads in phase (register 0), abundance 10 each
    positions = [(100, 10), (121, 10), (142, 10)]
    score, number, ratio = phase_metrics(positions, region_start=100, phase_length=21)
    assert number == 3
    assert ratio == 1.0
    assert score == round(1.0 * 3 * math.log(30), 3)

def test_phase_metrics_no_reads():
    assert phase_metrics([], 0, 21) == (0.0, 0, 0.0)

def test_hypergeometric_pvalue_bounds():
    p = hypergeometric_phase_pvalue(population=41, phased_positions=21, draws=5, observed=5)
    assert 0.0 <= p <= 1.0
    # observing 0 in phase is never significant
    assert hypergeometric_phase_pvalue(41, 21, 5, 0) == 1.0

def test_hypergeometric_degenerate_does_not_crash():
    # tiny cluster: population < phased_positions must return 1.0, not crash
    # (regression for math.comb ValueError on a negative argument)
    assert hypergeometric_phase_pvalue(population=9, phased_positions=27, draws=3, observed=2) == 1.0
    # more draws than the population is also degenerate
    assert hypergeometric_phase_pvalue(population=5, phased_positions=3, draws=10, observed=2) == 1.0
