"""5' first-base A fraction.

In plants the 5' nucleotide sorts a small RNA into an Argonaute: 5'U goes to
AGO1 (miRNA, tasiRNA/phasiRNA, ~21nt) and 5'A to AGO4 (heterochromatic siRNA,
24nt, the RdDM pathway). pct_1u already covers the U side; pct_5a is the A side.
Not to be confused with pct_10a, which is the TENTH base and a ping-pong
signature.
"""
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _aln(seq, pos=100, strand="+", nh=1):
    return Aln("c", pos, strand, seq, nh)


def test_pct_5a_counts_first_base_a():
    reads = [_aln("A" + "C" * 23), _aln("A" + "G" * 23), _aln("T" + "C" * 23),
             _aln("C" + "C" * 23)]
    r = build_record(("c", 100, 200), reads, 1_000_000, CFG)
    assert r.pct_5a == 0.5
    assert r.pct_1u == 0.25


def test_pct_5a_uses_the_oriented_read():
    # a minus-strand read is reverse-complemented before the base is read, so
    # this one starts with A in its own orientation
    minus = _aln("C" * 23 + "T", strand="-")
    r = build_record(("c", 100, 200), [minus], 1_000_000, CFG)
    assert r.pct_5a == 1.0


def test_pct_5a_is_weighted_by_nh():
    reads = [_aln("A" + "C" * 23, nh=4), _aln("T" + "C" * 23, nh=1)]
    r = build_record(("c", 100, 200), reads, 1_000_000, CFG)
    assert r.pct_5a == 0.2      # 0.25 / (0.25 + 1)
