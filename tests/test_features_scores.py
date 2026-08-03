"""build_record must turn the raw features into evidence scores."""
from features import build_record
from samio import Aln

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _reads(seq, n, pos=100, strand="+"):
    return [Aln("c", pos + i, strand, seq, 1) for i in range(n)]


def test_strong_5a_bias_scores_high():
    r = build_record(("c", 100, 300), _reads("A" + "C" * 23, 200),
                     1_000_000, CFG)
    assert r.s_5a > 20
    assert r.s_1u == 0.0


def test_shallow_locus_scores_low_despite_same_proportion():
    deep = build_record(("c", 100, 300), _reads("A" + "C" * 23, 200),
                        1_000_000, CFG)
    shallow = build_record(("c", 100, 300), _reads("A" + "C" * 23, 5),
                           1_000_000, CFG)
    assert shallow.s_5a < deep.s_5a


def test_single_strand_locus_scores_strand_bias():
    r = build_record(("c", 100, 300), _reads("A" + "C" * 23, 200),
                     1_000_000, CFG)
    assert r.s_fractop > 10


def test_phase_score_mirrors_the_pvalue():
    from evidence import score_from_p
    reads = [Aln("c", 100 + i * 21, "+", "T" + "A" * 20, 1) for i in range(11)]
    r = build_record(("c", 100, 400), reads, 1_000_000, CFG)
    assert r.s_phase == score_from_p(r.phase_pvalue)


def test_density_scored_only_when_genome_size_given():
    reads = _reads("A" + "C" * 23, 200)
    without = build_record(("c", 100, 300), reads, 1_000_000, CFG)
    with_genome = build_record(("c", 100, 300), reads, 1_000_000, CFG,
                               genome_bp=10_000_000)
    assert without.s_density == 0.0
    assert with_genome.s_density > 0
