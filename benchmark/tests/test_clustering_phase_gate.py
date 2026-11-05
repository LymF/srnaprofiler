"""Task 12: aggregate merge profile, strand_min_reads, nested windows, phase gate."""
import math
from collections import namedtuple
from types import SimpleNamespace

import pytest

from classify import class_size_windows
from clustering import merge_compatible
from decider import RuleScorer
from features import island_profiles

Aln = namedtuple("Aln", "pos seq nh strand")
CAln = namedtuple("CAln", "chrom pos seq nh strand")

PAD = 50


def _island(start, plus, minus, length=22):
    """One island at `start` carrying `plus`/`minus` reads (nh=1)."""
    alns = [Aln(start, "A" * length, 1, "+")] * plus
    alns += [Aln(start, "A" * length, 1, "-")] * minus
    return (start, start + length), alns


def _merge(specs, strand_min_reads, profile_min_reads=0, size_windows=None):
    """Merged loci for islands 100 bp apart (gap 78 < 2*PAD)."""
    islands, alns = [], []
    for k, (plus, minus) in enumerate(specs):
        isl, a = _island(k * 100, plus, minus)
        islands.append(isl)
        alns += a
    profiles = island_profiles(islands, alns, min_reads=profile_min_reads,
                               size_windows=size_windows)
    return merge_compatible(islands, profiles, PAD,
                            strand_min_reads=strand_min_reads,
                            profile_min_reads=profile_min_reads,
                            size_windows=size_windows)


# --- aggregate running profile ------------------------------------------------

def test_island_profiles_report_masses_for_aggregation():
    islands = [(0, 22), (100, 122)]
    _, a0 = _island(0, 3, 1)
    _, a1 = _island(100, 1, 0)
    prof = island_profiles(islands, a0 + a1, min_reads=2, size_windows=[(21, 22)])
    assert prof[0]["mass"] == 4 and prof[0]["plus"] == 3
    assert prof[0]["window_mass"] == [4]
    assert prof[0]["frac_top"] == 0.75 and prof[0]["window"] == 0
    # below min_reads: profile unknown, but the mass is still reported
    assert prof[1]["frac_top"] is None and prof[1]["window"] is None
    assert prof[1]["mass"] == 1 and prof[1]["plus"] == 1


def test_alternating_strands_merge_into_one_locus_when_strand_never_blocks():
    specs = [(30, 0), (0, 30), (30, 0), (0, 30)]
    assert len(_merge(specs, strand_min_reads=None)) == 1


def test_alternating_strands_split_when_both_sides_deep_enough():
    specs = [(30, 0), (0, 30), (30, 0), (0, 30)]
    assert len(_merge(specs, strand_min_reads=20)) == 4


def test_new_island_is_compared_with_aggregate_not_first_island():
    # first: 30+ (frac 1.0, '+'); second 10- is too thin to block and joins ->
    # aggregate 0.75 ('.'). Third is 30- (0.0, '-'): opposite to the first island
    # alone, but the aggregate is unstranded, so it merges.
    assert len(_merge([(30, 0), (0, 10), (0, 30)], strand_min_reads=20)) == 1


def test_aggregate_can_also_block_where_first_island_would_allow():
    # first 0.5 ('.'); second 60- allowed against it; aggregate 15/90 = 0.17 ('-').
    # Third 30+ (1.0): vs the first island it is unstranded and would merge, vs
    # the aggregate it is opposite and is blocked.
    out = _merge([(15, 15), (0, 60), (30, 0)], strand_min_reads=20)
    assert len(out) == 2


# --- strand_min_reads ---------------------------------------------------------

def test_strand_min_reads_thin_opposite_islands_merge():
    assert len(_merge([(5, 0), (0, 5)], strand_min_reads=20)) == 1


def test_strand_min_reads_deep_opposite_islands_do_not_merge():
    assert len(_merge([(50, 0), (0, 50)], strand_min_reads=20)) == 2


def test_strand_blocks_if_only_one_side_is_thin():
    assert len(_merge([(50, 0), (0, 5)], strand_min_reads=20)) == 1
    assert len(_merge([(5, 0), (0, 50)], strand_min_reads=20)) == 1


def test_island_below_profile_min_reads_never_blocks_and_adds_mass():
    # I1 (4 minus) is below profile_min_reads=5: unknown, joins I0 even with
    # strand_min_reads=0. Its mass still enters the aggregate (10+ / 4-, 0.714,
    # '.'), so I2 (1+/19-, 0.05, '-') merges; against I0 alone ('+') it would not.
    out = _merge([(10, 0), (0, 4), (1, 19)], strand_min_reads=0, profile_min_reads=5)
    assert len(out) == 1


# --- nested / duplicate size windows ------------------------------------------

def _cfg(**enabled):
    spec = lambda lo, hi, on=True: {"enabled": on, "size": {"lo": lo, "hi": hi}}
    return {"classes": {
        "phasiRNA21": spec(21, 22), "hcsiRNA": spec(23, 24),
        "siRNA": spec(21, 22),
        "phasiRNA24": spec(24, 24, enabled.get("p24", True))}}


def test_class_size_windows_drops_duplicate_and_nested_windows():
    assert class_size_windows(_cfg()) == [(21, 22), (23, 24)]


def test_enabling_nested_class_does_not_change_windows_or_merging():
    on, off = class_size_windows(_cfg()), class_size_windows(_cfg(p24=False))
    assert on == off
    # a 24-only island next to a 23/24 mix stays one locus either way
    islands = [(0, 24), (100, 124)]
    alns = [Aln(0, "A" * 24, 1, "+")] * 10
    alns += [Aln(100, "A" * 23, 1, "+")] * 5 + [Aln(100, "A" * 24, 1, "+")] * 5
    for w in (on, off):
        prof = island_profiles(islands, alns, size_windows=w)
        assert len(merge_compatible(islands, prof, PAD, size_windows=w)) == 1


def test_partially_overlapping_windows_are_kept():
    cfg = {"classes": {"a": {"size": {"lo": 20, "hi": 24}},
                       "b": {"size": {"lo": 23, "hi": 31}}}}
    assert class_size_windows(cfg) == [(20, 24), (23, 31)]


# --- phase_max_pvalue ---------------------------------------------------------

def _scorer(max_p):
    spec = {"size": {"lo": 24, "hi": 24}, "w_size": 1.0, "w_phase": 4.0,
            "phase_on": "cycles", "defining_evidence": ["phase"]}
    if max_p is not None:
        spec["phase_max_pvalue"] = max_p
    return RuleScorer({"classes": {"c": spec},
                       "decision": {"review_confidence": 0.5, "min_margin": 0.0,
                                    "min_reads": 0}}), spec


def _rec(p):
    return SimpleNamespace(dicer_call="24", phase_number=10, phase_score=50.0,
                           phase_pvalue=p, e_phase=0.8)


def test_phase_gate_zeroes_term_and_defining_strength():
    sc, spec = _scorer(0.01)
    _, terms = sc._score_class("c", spec, _rec(0.5))
    assert terms["phase"] == 0
    assert sc.defining_scores(spec, _rec(0.5)) == [0.0]


def test_phase_gate_passes_significant_record():
    sc, spec = _scorer(0.01)
    _, terms = sc._score_class("c", spec, _rec(0.001))
    assert terms["phase"] > 0
    assert sc.defining_scores(spec, _rec(0.001)) == [0.8]


@pytest.mark.parametrize("mode", ["cycles", "score"])
def test_phase_gate_applies_in_both_modes(mode):
    sc, spec = _scorer(0.01)
    spec["phase_on"] = mode
    assert sc._score_class("c", spec, _rec(0.5))[1]["phase"] == 0
    assert sc._score_class("c", spec, _rec(0.001))[1]["phase"] > 0


def test_phase_gate_nan_or_missing_pvalue_counts_as_zero():
    sc, spec = _scorer(0.01)
    assert sc._score_class("c", spec, _rec(float("nan")))[1]["phase"] == 0
    assert sc._score_class("c", spec, _rec(None))[1]["phase"] == 0
    rec = _rec(0.0)
    del rec.phase_pvalue
    assert sc.defining_scores(spec, rec) == [0.0]


def test_no_gate_when_key_absent():
    sc, spec = _scorer(None)
    assert sc._score_class("c", spec, _rec(0.5))[1]["phase"] > 0
    assert sc._score_class("c", spec, _rec(float("nan")))[1]["phase"] > 0
    assert sc.defining_scores(spec, _rec(0.5)) == [0.8]


def test_phase_gate_does_not_touch_other_defining_evidence():
    sc, spec = _scorer(0.01)
    spec["defining_evidence"] = [["phase", "5a"]]
    rec = _rec(0.5)
    rec.e_5a = 0.3
    assert sc.defining_scores(spec, rec) == [0.3]


# --- strand call (ShortStack strand_cutoff) -----------------------------------

def test_opposite_stranded_islands_do_not_merge():
    # frac_top 0.9 ('+') vs 0.1 ('-') with 30 reads each
    assert len(_merge([(27, 3), (3, 27)], strand_min_reads=20)) == 2


def test_stranded_next_to_unstranded_merges():
    assert len(_merge([(27, 3), (15, 15)], strand_min_reads=20)) == 1
    assert len(_merge([(15, 15), (3, 27)], strand_min_reads=20)) == 1


def test_opposite_stranded_but_thin_side_merges():
    assert len(_merge([(27, 3), (1, 9)], strand_min_reads=20)) == 1


def test_strand_cutoff_is_a_parameter_and_a_boundary_is_stranded():
    # 0.8 and 0.2 sit exactly on the cutoff: stranded (>= c, <= 1 - c)
    assert len(_merge([(24, 6), (6, 24)], strand_min_reads=20)) == 2
    islands, alns = [], []
    for k, (p, m) in enumerate([(27, 3), (3, 27)]):
        isl, a = _island(k * 100, p, m)
        islands.append(isl)
        alns += a
    prof = island_profiles(islands, alns)
    # 0.9 / 0.1 are not stranded under a 0.95 cutoff
    assert len(merge_compatible(islands, prof, PAD, strand_min_reads=20,
                                strand_cutoff=0.95)) == 1


def test_process_chromosome_reads_cutoff_from_config(monkeypatch):
    import features
    seen = {}

    def fake_merge(islands, profiles, pad, **kw):
        seen.update(kw)
        return [(0, 100)]
    alns = [CAln("chr1", a.pos, a.seq, a.nh, a.strand) for a in _island(0, 3, 0)[1]]
    monkeypatch.setattr(features, "read_region", lambda *a, **k: iter(alns))
    monkeypatch.setattr(features, "merge_compatible", fake_merge)
    monkeypatch.setattr(features, "build_record",
                        lambda region, *a, **k: SimpleNamespace(
                            start=region[1], end=region[2], reads=3.0))
    cfg = {"strand_cutoff": 0.77, "strand_min_reads": 20}
    features.process_chromosome("x.bam", "chr1", 1, 50, 1e6, cfg, chrom_len=1000)
    assert seen["strand_cutoff"] == 0.77


# --- gated classes leave the softmax ------------------------------------------

def _p24_scorer():
    cfg = {"classes": {
        "a": {"size": {"lo": 24, "hi": 24}, "w_size": 1.0},
        "b": {"size": {"lo": 24, "hi": 24}, "w_size": 0.5},
        "p24": {"size": {"lo": 24, "hi": 24}, "w_size": 1.0, "w_phase": 4.0,
                "phase_on": "cycles", "phase_max_pvalue": 0.01,
                "defining_evidence": ["phase"]}},
        "decision": {"review_confidence": 0.5, "min_margin": 0.0, "min_reads": 0}}
    return RuleScorer(cfg)


def _drec(p, reads=100):
    return SimpleNamespace(dicer_call="24", phase_number=10, phase_score=50.0,
                           phase_pvalue=p, e_phase=0.8, reads=reads,
                           size_counts={}, distinct_seqs=0)


def test_failing_gate_class_leaves_the_softmax():
    rec = _p24_scorer().decide(_drec(0.5))
    assert rec.cls == "a"
    expected = math.exp(1.0) / (math.exp(1.0) + math.exp(0.5))
    assert rec.confidence == pytest.approx(expected, abs=1e-3)
    assert rec.class_scores["p24"] is None
    assert rec.class_scores["a"] == 1.0


def test_passing_gate_class_competes():
    rec = _p24_scorer().decide(_drec(0.001))
    assert rec.cls == "p24"
    assert rec.class_scores["p24"] > 1.0


def test_all_gated_classes_failing_falls_back_to_scoring_all():
    cfg = {"classes": {"p": {"size": {"lo": 24, "hi": 24}, "w_size": 1.0,
                             "w_phase": 1.0, "phase_max_pvalue": 0.01}},
           "decision": {"review_confidence": 0.5, "min_margin": 0.0,
                        "min_reads": 0}}
    rec = RuleScorer(cfg).decide(_drec(0.5))
    assert rec.cls == "p"
    assert rec.class_scores["p"] == 1.0


def test_tsv_writer_handles_excluded_class(tmp_path):
    from reporter import write_tsv
    from records import ClusterRecord
    r = ClusterRecord(chrom="c", start=0, end=10)
    r.class_scores = {"x": None}
    write_tsv([r], str(tmp_path / "o.tsv"))


# --- proTRAC piRNA-cluster criterion ------------------------------------------

def _pi_scorer(**keys):
    spec = {"size": {"lo": 24, "hi": 31}, "w_size": 1.0}
    spec.update(keys)
    other = {"size": {"lo": 21, "hi": 22}, "w_size": 0.5}
    return RuleScorer({"classes": {"piRNA": spec, "siRNA": other},
                       "decision": {"review_confidence": 0.0, "min_margin": 0.0,
                                    "min_reads": 0}})


def _crec(length, dens, dicer="27"):
    return SimpleNamespace(dicer_call=dicer, reads=100, start=0, end=length,
                           density_pvalue=dens, size_counts={}, distinct_seqs=0,
                           phase_pvalue=1.0, phase_number=0)


KEYS = dict(cluster_min_length=1000, cluster_max_density_pvalue=0.01)


def test_short_pirna_locus_is_downgraded():
    rec = _pi_scorer(**KEYS).decide(_crec(800, 0.001))
    assert rec.cls == "ambiguous" and rec.review is True and rec.score is None
    assert rec.evidence.startswith("nearest=piRNA; below_cluster_criterion(length); ")


def test_long_dense_pirna_locus_stays():
    rec = _pi_scorer(**KEYS).decide(_crec(1500, 0.001))
    assert rec.cls == "piRNA"


def test_long_sparse_pirna_locus_is_downgraded():
    rec = _pi_scorer(**KEYS).decide(_crec(1500, 0.5))
    assert rec.cls == "ambiguous" and rec.review is True
    assert "below_cluster_criterion(density)" in rec.evidence


def test_both_failures_are_named():
    rec = _pi_scorer(**KEYS).decide(_crec(800, 0.5))
    assert "below_cluster_criterion(length,density)" in rec.evidence


def test_class_without_keys_is_unaffected():
    rec = _pi_scorer().decide(_crec(800, 0.5))
    assert rec.cls == "piRNA"


def test_missing_density_applies_only_length():
    sc = _pi_scorer(**KEYS)
    assert sc.decide(_crec(1500, None)).cls == "piRNA"
    rec = _crec(1500, None)
    del rec.density_pvalue
    assert sc.decide(rec).cls == "piRNA"
    assert sc.decide(_crec(800, None)).cls == "ambiguous"


def test_nan_density_counts_as_missing():
    assert _pi_scorer(**KEYS).decide(_crec(1500, float("nan"))).cls == "piRNA"


def test_criterion_never_changes_the_winner_only_downgrades_it():
    rec = _pi_scorer(**KEYS).decide(_crec(800, 0.5, dicer="21"))
    assert rec.cls == "siRNA"


# --- density_pvalue -----------------------------------------------------------

def test_density_pvalue_hand_checked():
    from features import density_pvalue
    # chromosome: 1000 reads over 100000 bp; locus 200 bp -> lambda = 2.
    # P(X >= 3) = 1 - e^-2 (1 + 2 + 2) = 1 - 5 e^-2
    assert density_pvalue(3, 1000, 100000, 200) == pytest.approx(
        1 - 5 * math.exp(-2), abs=1e-9)
    # integer k agrees with the Poisson tail; non-integer is the gamma extension
    assert density_pvalue(0, 1000, 100000, 200) == 1.0
    assert density_pvalue(2.5, 1000, 100000, 200) < density_pvalue(2, 1000, 100000, 200)
    assert density_pvalue(3, 1000, 0, 200) is None


def test_process_chromosome_stores_density_pvalue(monkeypatch):
    import features
    alns = [CAln("chr1", a.pos, a.seq, a.nh, a.strand) for a in _island(0, 3, 0)[1]]
    monkeypatch.setattr(features, "read_region", lambda *a, **k: iter(alns))
    monkeypatch.setattr(features, "merge_compatible", lambda *a, **k: [(0, 200)])
    monkeypatch.setattr(features, "build_record",
                        lambda region, *a, **k: SimpleNamespace(
                            start=region[1], end=region[2], reads=3.0))
    cfg = {"strand_cutoff": 0.8}
    # chrom mass 3 over 300 bp -> 0.01/bp; locus 200 bp -> lambda 2; k = 3
    recs = features.process_chromosome("x.bam", "chr1", 1, 50, 1e6, cfg,
                                       chrom_len=300)
    assert recs[0].density_pvalue == pytest.approx(1 - 5 * math.exp(-2), abs=1e-9)
