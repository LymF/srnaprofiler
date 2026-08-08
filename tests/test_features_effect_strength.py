"""Each defining evidence gets a strength alongside its p-value.

They are not TSV columns: the decider reads them to build `score`.
"""
from itertools import product

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
