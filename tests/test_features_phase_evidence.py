"""e_phase measures how far the phasing runs, not how surprising it is.

Measured on the 8 experimental TAS (12 loci) and 68 loci on annotated PHAS21,
against controls matched by distinct_seqs: the count of cycles in register is
the best discriminator available (AUC 0.826 / 0.769), the current p-value is at
chance (0.560 / 0.481) and the abundance ratio is worse than a coin (0.415).
Dividing the count by the locus's own occupied cycles destroys it (0.647),
because the positives simply span more cycles than their controls.
"""
from evidence import effect_strength
from features import build_record
from samio import Aln
from signals import PHASE_WINDOW_CYCLES

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}

_TAILS = ["ACGT", "TGCA", "GTAC", "CATG", "AGCT", "TCGA", "GATC", "CTAG",
          "AACC", "GGTT", "CCAA", "TTGG", "ACAC", "GTGT", "CACA", "TGTG"]


def _seqs(n, length=21):
    return [(_TAILS[i % len(_TAILS)] * ((length // 4) + 1))[:length]
            for i in range(n)]


def _in_register(n_cycles, phase=21, start=0):
    """One 21-nt read per phase cycle, all in the same register."""
    return [Aln("c", start + i * phase, "+", seq, 1)
            for i, seq in enumerate(_seqs(n_cycles, phase))]


def test_a_full_window_in_register_is_the_maximum():
    r = build_record(("c", 0, 400), _in_register(PHASE_WINDOW_CYCLES),
                     1_000_000, CFG)
    assert r.dicer_call == "21"
    assert r.e_phase == effect_strength(PHASE_WINDOW_CYCLES,
                                        PHASE_WINDOW_CYCLES, 1 / 21)
    assert r.e_phase > 0.6


def test_more_cycles_in_register_is_more_evidence():
    """The whole point: the statistic has to grow with the run's length."""
    short = build_record(("c", 0, 400), _in_register(3), 1_000_000, CFG)
    long = build_record(("c", 0, 400), _in_register(9), 1_000_000, CFG)
    assert long.e_phase > short.e_phase


def test_the_denominator_is_the_window_not_the_locus():
    """A 3-cycle locus must not score like a 3-cycle run inside a long one.

    Normalising by the locus's own cycles was measured and rejected: it drops
    the AUC from 0.826 to 0.647, because a TAS spans 10 cycles and its matched
    control spans 5.
    """
    r = build_record(("c", 0, 400), _in_register(3), 1_000_000, CFG)
    assert r.e_phase == effect_strength(3, PHASE_WINDOW_CYCLES, 1 / 21)
    assert r.e_phase < 0.3


def test_abundance_alone_does_not_make_evidence():
    """One hot position with 40 reads is stacking, not processivity."""
    stacked = [Aln("c", 100, "+", seq, 1) for seq in _seqs(40)]
    r = build_record(("c", 0, 400), stacked, 1_000_000, CFG)
    assert r.e_phase < 0.3


def test_an_empty_locus_has_no_phase_evidence():
    r = build_record(("c", 0, 400), [], 1_000_000, CFG)
    assert r.e_phase == 0.0


def test_a_locus_with_no_dicer_call_has_no_phase_evidence():
    """Sizes all over the window: there is no phase length to speak of."""
    seqs = _seqs(15, 40)
    # Mix sizes: 5 short, 5 in-window, 5 long to ensure >20% outside [18, 35]
    mixed = (
        [Aln("c", i * 10, "+", seq[:16], 1) for i, seq in enumerate(seqs[:5])] +
        [Aln("c", i * 10 + 50, "+", seq[:27], 1) for i, seq in enumerate(seqs[5:10])] +
        [Aln("c", i * 10 + 100, "+", seq[:37], 1) for i, seq in enumerate(seqs[10:15])]
    )
    r = build_record(("c", 0, 900), mixed, 1_000_000, CFG)
    assert r.dicer_call == "N"
    assert r.e_phase == 0.0
