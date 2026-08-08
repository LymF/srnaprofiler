"""Each defining evidence gets a strength alongside its p-value.

They are not TSV columns: the decider reads them to build `score`.
"""
from itertools import product

from evidence import effect_strength
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}

_TAILS = ["".join(p) for p in product("ACGT", repeat=4)]


def _reads(first_base, n, start=100, strand="+"):
    """n reads of 24 nt, all starting with `first_base`, each sequence unique.

    Unique because n is distinct_seqs: a stack of identical reads is one
    observation, which is the pseudoreplication guard the tool applies on purpose.
    """
    out = []
    for i in range(n):
        seq = (first_base + _TAILS[i % len(_TAILS)] + "ACGT" * 5)[:24]
        out.append(Aln("c", start + i * 3, strand, seq, 1))
    return out


def test_a_strong_1u_locus_gets_a_high_strength():
    r = build_record(("c", 0, 400), _reads("T", 60), 1_000_000, CFG)
    assert r.e_1u > 0.5


def test_a_locus_without_1u_gets_zero():
    r = build_record(("c", 0, 400), _reads("G", 60), 1_000_000, CFG)
    assert r.e_1u == 0.0


def test_strength_does_not_grow_with_the_locus():
    """The p-value score does; that is the bug this replaces."""
    small = build_record(("c", 0, 400), _reads("T", 40), 1_000_000, CFG)
    big = build_record(("c", 0, 4000), _reads("T", 400), 1_000_000, CFG)
    assert abs(big.e_1u - small.e_1u) < 0.25
    assert big.s_1u > small.s_1u * 3      # the p-value score, for contrast


def test_bidirectional_strength_needs_both_strands():
    uni = build_record(("c", 0, 400), _reads("T", 60), 1_000_000, CFG)
    dual = build_record(("c", 0, 400),
                        _reads("T", 30)
                        + _reads("T", 30, start=200, strand="-"),
                        1_000_000, CFG)
    assert dual.e_fractop_bid > uni.e_fractop_bid
    assert uni.e_fractop_bid == 0.0


def test_an_empty_locus_has_zero_strengths():
    r = build_record(("c", 0, 400), [], 1_000_000, CFG)
    for attr in ("e_1u", "e_5a", "e_fractop_bid", "e_pingpong", "e_phase"):
        assert getattr(r, attr) == 0.0


def _distinct_seqs(n, length=21):
    """n distinct sequences of `length` nt, built from the 4-nt tail pool."""
    return [(_TAILS[i % len(_TAILS)] * ((length // 4) + 1))[:length] for i in range(n)]


def _at(position, strand, n, length=21):
    """n reads of `length` nt at a fixed genomic position, distinct sequences."""
    return [Aln("c", position, strand, seq, 1) for seq in _distinct_seqs(n, length)]


def test_a_strong_10nt_overlap_gets_a_high_pingpong_strength():
    # minus 5' end = pos + length - 1, so overlap d needs minus.pos = plus.pos + d - length.
    # 10 plus reads at 100 pair with 10 minus reads at 89 (5'-5' overlap = 10),
    # plus a few off-peak pairs so the background is not all zero.
    reads = (_at(100, "+", 10)
             + _at(89, "-", 10)     # overlap 10, the peak
             + _at(84, "-", 3)      # overlap 5, background
             + _at(94, "-", 3))     # overlap 15, background
    r = build_record(("c", 0, 400), reads, 1_000_000, CFG)
    assert r.e_pingpong > 0.3


def test_a_clean_phase_register_gets_a_high_phase_strength():
    # 11 reads of 21 nt, one per cycle, in perfect register -> dicer_call "21".
    from signals import PHASE_WINDOW_CYCLES

    reads = [Aln("c", i * 21, "+", seq, 1)
             for i, seq in enumerate(_distinct_seqs(11, 21))]
    r = build_record(("c", 0, 400), reads, 1_000_000, CFG)
    assert r.dicer_call == "21"
    # A full window in register is as much phase evidence as there is.
    assert r.e_phase == effect_strength(PHASE_WINDOW_CYCLES,
                                        PHASE_WINDOW_CYCLES, 1 / 21)


def test_a_pure_10nt_overlap_with_no_background_gets_a_positive_strength():
    """Only offset 10 is populated -- the exact case the zero-background guard broke."""
    reads = _at(100, "+", 10) + _at(89, "-", 10)   # overlap 10 only, no other offsets
    r = build_record(("c", 0, 400), reads, 1_000_000, CFG)
    assert r.e_pingpong > 0
