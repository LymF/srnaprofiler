"""build_record must turn the raw features into evidence scores."""
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _reads(seq, n, pos=100, strand="+"):
    """n copies of the SAME sequence (n_distinct_seqs == 1): stacked reads."""
    return [Aln("c", pos + i, strand, seq, 1) for i in range(n)]


def _distinct_reads(n, pos=100, strand="+", base="A"):
    """n DISTINCT sequences sharing a 5' base: independent trials, not a stack."""
    return [Aln("c", pos + i, strand, base + f"{i:03d}C" * 6, 1) for i in range(n)]


def test_strong_5a_bias_scores_high():
    r = build_record(("c", 100, 300), _distinct_reads(200, base="A"),
                     1_000_000, CFG)
    assert r.s_5a > 20
    assert r.s_1u == 0.0


def test_shallow_locus_scores_low_despite_same_proportion():
    deep = build_record(("c", 100, 300), _distinct_reads(200, base="A"),
                        1_000_000, CFG)
    shallow = build_record(("c", 100, 300), _distinct_reads(5, base="A"),
                           1_000_000, CFG)
    assert shallow.s_5a < deep.s_5a


def test_single_strand_locus_scores_strand_bias():
    r = build_record(("c", 100, 300), _distinct_reads(200, base="A"),
                     1_000_000, CFG)
    assert r.s_fractop > 10


def test_stacked_identical_reads_do_not_inflate_the_score():
    # pseudoreplication guard: n must be distinct_seqs, not the 1/NH read count,
    # so a pile of IDENTICAL reads (one distinct sequence) cannot fake significance.
    r = build_record(("c", 100, 300), _reads("A" + "C" * 23, 200), 1_000_000, CFG)
    assert r.distinct_seqs == 1
    assert r.s_5a < 1.0


def test_dual_strand_locus_scores_high_on_fractop_bal():
    # frac_top ~0.5 (dual-strand, siRNA-defining) must score HIGH on s_fractop_bal
    # and LOW on s_fractop -- the opposite pattern, since the two null tests differ.
    top = _distinct_reads(100, base="A", strand="+", pos=100)
    bottom = _distinct_reads(100, base="A", strand="-", pos=200)
    r = build_record(("c", 100, 400), top + bottom, 1_000_000, CFG)
    assert r.frac_top == 0.5
    assert r.s_fractop_bal > 5
    assert r.s_fractop < 0.5


def test_uni_strand_locus_scores_low_on_fractop_bal():
    r = build_record(("c", 100, 300), _distinct_reads(200, base="A"),
                     1_000_000, CFG)
    assert r.s_fractop_bal < 0.5


def test_phase_score_mirrors_the_pvalue():
    from evidence import score_from_p
    reads = [Aln("c", 100 + i * 21, "+", "T" + "A" * 20, 1) for i in range(11)]
    r = build_record(("c", 100, 400), reads, 1_000_000, CFG)
    assert r.s_phase == score_from_p(r.phase_pvalue)


def test_density_is_no_longer_reported():
    """It restated the read count (r=0.72 with log reads) and was high by
    construction, since a locus only exists after clearing a depth cut."""
    r = build_record(("c", 100, 300), _reads("A" + "C" * 23, 200), 1_000_000, CFG)
    assert not hasattr(r, "s_density")
