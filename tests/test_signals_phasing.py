import math
from signals import phase_metrics, hypergeometric_phase_pvalue

def test_phase_metrics_perfectly_phased():
    # três reads em fase (registro 0), abundância 10 cada
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
    # observar 0 em fase nunca é significativo
    assert hypergeometric_phase_pvalue(41, 21, 5, 0) == 1.0
