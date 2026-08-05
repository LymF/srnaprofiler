"""The Dicer duplex peak sits at the locus's own size minus 2, not at a fixed 19.

Two strands of one duplex leave their 5' ends L-2 apart (2-nt 3' overhang). The
old column read offset 19 for every locus, which is only right at 21 nt and was
outside the spectrum window for anything above 22.
"""
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _duplex(length, n, start=100, step=40):
    """n duplexes of `length`, each a + read and its - partner (5' ends L-2 apart)."""
    alns = []
    for i in range(n):
        p = start + i * step
        alns.append(Aln("c", p, "+", "A" * length, 1))
        # the - read's 5' end is its last base: place it so the overlap is L-2
        alns.append(Aln("c", p - 2, "-", "A" * length, 1))
    return alns


def test_duplex_detected_at_21nt():
    r = build_record(("c", 0, 3000), _duplex(21, 20), 1_000_000, CFG)
    assert r.dicer_call == "21"
    assert r.overlap_dup_frac > 0.5


def test_duplex_detected_at_27nt():
    # this is the case the fixed offset 19 could not see at all
    r = build_record(("c", 0, 3000), _duplex(27, 20), 1_000_000, CFG)
    assert r.dicer_call == "27"
    assert r.overlap_dup_frac > 0.5


def test_no_duplex_geometry_gives_no_signal():
    # all reads on one strand: no cross-strand pairs at any offset
    alns = [Aln("c", 100 + i * 40, "+", "A" * 21, 1) for i in range(20)]
    r = build_record(("c", 0, 3000), alns, 1_000_000, CFG)
    assert r.overlap_dup_frac == 0.0


def test_unknown_size_gives_no_duplex_call():
    r = build_record(("c", 0, 3000), [], 1_000_000, CFG)
    assert r.dicer_call == "N"
    assert r.overlap_dup_frac == 0.0


def test_pingpong_fraction_is_unchanged_by_the_wider_window():
    # 5'-5' overlap of 10: ping-pong geometry, must keep its own 1..20 denominator
    alns = []
    for i in range(20):
        p = 100 + i * 60
        alns.append(Aln("c", p, "+", "A" * 26, 1))
        alns.append(Aln("c", p - 16, "-", "A" * 26, 1))   # five_minus = p + 9
    r = build_record(("c", 0, 3000), alns, 1_000_000, CFG)
    assert r.pingpong_frac > 0.5
