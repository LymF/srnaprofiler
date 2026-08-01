from samio import Aln
from class_signals import base_composition


def test_base_composition_weights_and_orients():
    alns = [
        Aln("chr", 0, "+", "TACGT", 1),       # len 5, + -> first base T, w=1
        Aln("chr", 0, "-", "AACGT", 2),        # revcomp = ACGTT -> first base A, w=0.5
    ]
    comp = base_composition(alns, minlen=5, maxlen=5)
    assert comp[5]['+']['T'] == 1.0
    assert comp[5]['-']['A'] == 0.5
    assert comp[5]['+']['A'] == 0.0


def test_base_composition_size_filter():
    alns = [Aln("chr", 0, "+", "AAA", 1)]      # len 3, outside 5..35
    comp = base_composition(alns, minlen=5, maxlen=35)
    assert all(comp[s]['+']['A'] == 0.0 for s in comp)
