import pytest

from benchmark.manifest import Config
from benchmark.metrics_class import classification_metrics, match_loci, read_loci
from benchmark.truth import Region

CFG = Config(tiny_bp=[50, 100], low_reads=[20, 50, 100], min_overlap_frac=0.5)


def R(start, end, cls, chrom="c1", name="r"):
    return Region(chrom, start, end, cls, "t", name)


def L(start, end, cls, conf=0.9, reads=10.0, chrom="c1"):
    return {"chrom": chrom, "start": start, "end": end, "cls": cls,
            "confidence": conf, "reads": reads}


def test_read_loci(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text("chrom\tstart\tend\treads\tcls\tconfidence\tname\n"
                 "c1\t0\t10\t5.5\tpiRNA\t0.8\tx\n"
                 "c1\t10\t20\t1\tambiguous\tNA\ty\n")
    loci = read_loci(p)
    assert loci[0]["start"] == 0 and loci[0]["reads"] == 5.5
    assert loci[0]["confidence"] == 0.8 and loci[0]["name"] == "x"
    assert loci[1]["confidence"] is None and loci[1]["end"] == 20


def test_read_loci_read_span(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text("chrom\tstart\tend\treads\tread_span\tcls\n"
                 "c1\t0\t10\t5\t8\tpiRNA\nc1\t10\t20\t1\tNA\tpiRNA\n")
    loci = read_loci(p)
    assert loci[0]["read_span"] == 8 and loci[1]["read_span"] is None
    p.write_text("chrom\tstart\tend\treads\tcls\nc1\t0\t10\t5\tpiRNA\n")
    assert read_loci(p)[0]["read_span"] is None


def test_match_best_overlap_and_threshold():
    regions = [R(0, 100, "piRNA", name="a"), R(50, 400, "siRNA", name="b")]
    out = match_loci([L(40, 200, "piRNA"), L(1000, 1100, "piRNA"),
                      L(0, 10, "piRNA", chrom="other")], regions, 0.5)
    assert out[0][1].name == "b"  # overlap 150 vs 60
    assert out[1][1] is None
    assert out[2][1] is None


def test_match_region_fraction_criterion():
    # locus much larger than region: overlap/region_len = 1.0
    out = match_loci([L(0, 1000, "piRNA")], [R(100, 150, "piRNA")], 0.5)
    assert out[0][1] is not None


def test_precision_recall_per_class():
    regions = [R(0, 100, "piRNA"), R(200, 300, "piRNA"), R(400, 500, "siRNA")]
    loci = [L(0, 100, "piRNA"), L(400, 500, "piRNA"), L(200, 300, "siRNA")]
    m = classification_metrics(loci, regions, CFG)["per_class"]
    assert m["piRNA"]["precision"] == pytest.approx(0.5)
    assert m["piRNA"]["fdr"] == pytest.approx(0.5)
    assert m["piRNA"]["recall"] == pytest.approx(0.5)
    assert m["siRNA"]["precision"] == 0.0
    assert m["siRNA"]["recall"] == 0.0


def test_ambiguous_is_abstention_not_error():
    regions = [R(0, 100, "piRNA")]
    loci = [L(0, 100, "piRNA"), L(0, 100, "ambiguous", conf=None)]
    res = classification_metrics(loci, regions, CFG)
    assert res["per_class"]["piRNA"]["precision"] == 1.0
    assert "ambiguous" not in res["per_class"]
    assert res["ambiguous_coverage"] == pytest.approx(0.5)


def test_unannotated_loci_reported_as_unknown_not_fp():
    regions = [R(0, 100, "piRNA")]
    loci = [L(0, 100, "piRNA"), L(1000, 1100, "piRNA", reads=20.0),
            L(2000, 2300, "ambiguous", conf=None, reads=40.0)]
    res = classification_metrics(loci, regions, CFG)
    assert res["per_class"]["piRNA"]["precision"] == 1.0
    assert res["unknown"]["count"] == 2
    assert res["unknown"]["length"]["p50"] == pytest.approx(200.0)
    assert res["unknown"]["reads"]["p50"] == pytest.approx(30.0)
    assert set(res["unknown"]["reads"]) == {"p5", "p25", "p50", "p75", "p95"}


def test_ece_perfect_calibration_is_zero():
    regions = [R(i * 100, i * 100 + 50, "piRNA") for i in range(10)]
    # confidence 0.5, half correct
    loci = [L(i * 100, i * 100 + 50, "piRNA" if i < 5 else "siRNA", conf=0.5)
            for i in range(10)]
    assert classification_metrics(loci, regions, CFG)["ece"] == pytest.approx(0.0)


def test_ece_overconfident():
    regions = [R(0, 50, "piRNA"), R(100, 150, "piRNA")]
    loci = [L(0, 50, "piRNA", conf=1.0), L(100, 150, "siRNA", conf=1.0)]
    assert classification_metrics(loci, regions, CFG)["ece"] == pytest.approx(0.5)


def test_empty_loci_returns_zero_metrics():
    res = classification_metrics([], [R(0, 100, "piRNA")], CFG)
    assert res["ambiguous_coverage"] == 0.0
    assert res["unknown"]["count"] == 0
    assert res["unknown"]["length"]["p50"] == 0.0
    assert res["ece"] == 0.0
    assert res["per_class"]["piRNA"]["precision"] is None
    assert res["per_class"]["piRNA"]["recall"] == 0.0


def test_match_below_threshold_on_both_fractions_is_unmatched():
    # overlap 10: 10/100 = 0.1 of locus, 10/910 ~ 0.011 of region
    out = match_loci([L(0, 100, "piRNA")], [R(90, 1000, "piRNA")], 0.5)
    assert out[0][1] is None


def test_match_exactly_at_threshold_matches():
    # overlap 50: 50/100 == 0.5 of locus; region fraction 0.05
    out = match_loci([L(0, 100, "piRNA")], [R(50, 1050, "piRNA")], 0.5)
    assert out[0][1] is not None


def test_match_long_region_starting_far_before_locus():
    regions = [R(0, 10000, "piRNA", name="long"), R(100, 150, "siRNA", name="s1"),
               R(5000, 5050, "siRNA", name="s2"), R(9500, 9550, "siRNA", name="s3")]
    out = match_loci([L(9000, 9100, "piRNA")], regions, 0.5)
    assert out[0][1].name == "long"


def test_ece_multiple_bins_different_sizes():
    regions = [R(i * 100, i * 100 + 50, "piRNA") for i in range(4)]
    loci = [L(0, 50, "piRNA", conf=0.9), L(100, 150, "piRNA", conf=0.9),
            L(200, 250, "siRNA", conf=0.9), L(300, 350, "siRNA", conf=0.2)]
    expected = 0.75 * abs(2 / 3 - 0.9) + 0.25 * abs(0.0 - 0.2)
    assert classification_metrics(loci, regions, CFG)["ece"] == pytest.approx(expected)


def test_ece_excludes_ambiguous_unmatched_and_missing_confidence():
    regions = [R(0, 50, "piRNA"), R(100, 150, "piRNA")]
    loci = [L(0, 50, "piRNA", conf=1.0),
            L(100, 150, "ambiguous", conf=0.1),
            L(1000, 1050, "piRNA", conf=0.1),
            L(100, 150, "siRNA", conf=None)]
    assert classification_metrics(loci, regions, CFG)["ece"] == pytest.approx(0.0)


def test_per_class_counts_called_matched_unknown():
    regions = [R(0, 100, "piRNA"), R(400, 500, "siRNA")]
    loci = [L(0, 100, "piRNA"),                      # tp
            L(400, 500, "piRNA"),                    # fp (matches an siRNA region)
            L(1000, 1100, "piRNA"),                  # called, no match -> unknown
            L(2000, 2100, "ambiguous", conf=None)]   # abstention: not called
    pc = classification_metrics(loci, regions, CFG)["per_class"]
    assert pc["piRNA"]["n_called"] == 3
    assert pc["piRNA"]["n_matched"] == 2 == pc["piRNA"]["tp"] + pc["piRNA"]["fp"]
    assert pc["piRNA"]["n_unknown"] == 1
    assert pc["piRNA"]["precision"] == pytest.approx(0.5)  # unknown loci are not in the denominator
    assert pc["siRNA"]["n_called"] == 0 and pc["siRNA"]["n_unknown"] == 0
    assert pc["piRNA"]["n_truth"] == 1 and pc["siRNA"]["n_truth"] == 1


def test_recall_expressed_uses_only_expressed_truth():
    regions = [R(0, 100, "piRNA", name="a"), R(200, 300, "piRNA", name="b"),
               R(400, 500, "piRNA", name="c")]
    loci = [L(0, 100, "piRNA"), L(400, 500, "piRNA")]
    expressed = {id(regions[0]), id(regions[1])}  # a, b expressed; c is silent
    pc = classification_metrics(loci, regions, CFG, expressed=expressed)["per_class"]["piRNA"]
    assert pc["recall"] == pytest.approx(2 / 3)
    assert pc["n_truth"] == 3 and pc["n_truth_expressed"] == 2
    assert pc["recall_expressed"] == pytest.approx(0.5)  # only a of {a, b} is hit


def test_recall_expressed_none_without_expressed_regions_or_bam():
    regions = [R(0, 100, "piRNA")]
    loci = [L(0, 100, "piRNA")]
    pc = classification_metrics(loci, regions, CFG)["per_class"]["piRNA"]
    assert pc["recall_expressed"] is None and pc["n_truth_expressed"] is None
    pc = classification_metrics(loci, regions, CFG, expressed=set())["per_class"]["piRNA"]
    assert pc["recall_expressed"] is None and pc["n_truth_expressed"] == 0


def test_match_tie_prefers_smaller_region_then_input_order():
    # equal overlap (100 bp) with a big region listed first and a small one second
    regions = [R(0, 1000, "siRNA", name="big"), R(0, 100, "miRNA", name="small")]
    out = match_loci([L(0, 100, "miRNA")], regions, 0.5)
    assert out[0][1].name == "small"
    regions = [R(0, 100, "siRNA", name="first"), R(0, 100, "miRNA", name="second")]
    out = match_loci([L(0, 100, "miRNA")], regions, 0.5)
    assert out[0][1].name == "first"
    # the order of the manifest does not change the winner among different sizes
    out = match_loci([L(0, 100, "miRNA")], list(reversed([R(0, 1000, "siRNA", name="big"),
                                                          R(0, 100, "miRNA", name="small")])), 0.5)
    assert out[0][1].name == "small"


def _pc(loci, regions, cls="piRNA", **kw):
    return classification_metrics(loci, regions, CFG, **kw)["per_class"][cls]


def test_recall_one_locus_covering_three_regions_recalls_all():
    regions = [R(0, 100, "piRNA", name="a"), R(110, 210, "piRNA", name="b"),
               R(220, 320, "piRNA", name="c")]
    pc = _pc([L(0, 320, "piRNA")], regions)
    assert pc["recall"] == pytest.approx(1.0)  # each region fully inside the locus
    assert pc["precision"] == pytest.approx(1.0)
    assert pc["tp"] == 1 and pc["fp"] == 0 and pc["n_called"] == 1  # per-locus counts unchanged


def test_recall_ignores_loci_of_other_class():
    regions = [R(0, 100, "piRNA")]
    m = classification_metrics([L(0, 100, "siRNA")], regions, CFG)["per_class"]
    assert m["piRNA"]["recall"] == 0.0


def test_recall_ignores_ambiguous_loci():
    regions = [R(0, 100, "piRNA")]
    assert _pc([L(0, 100, "ambiguous", conf=None)], regions)["recall"] == 0.0


def test_recall_overlap_criterion_both_directions():
    # small region inside a big locus: overlap/region_len = 1.0
    assert _pc([L(0, 1000, "piRNA")], [R(100, 150, "piRNA")])["recall"] == 1.0
    # small locus inside a big region: overlap/locus_len = 1.0
    assert _pc([L(100, 150, "piRNA")], [R(0, 1000, "piRNA")])["recall"] == 1.0
    # exactly 0.5 of the locus counts; just below on both fractions does not
    assert _pc([L(0, 100, "piRNA")], [R(50, 1000, "piRNA")])["recall"] == 1.0
    assert _pc([L(0, 100, "piRNA")], [R(51, 1000, "piRNA")])["recall"] == 0.0


def test_recall_region_lost_to_other_class_best_match_still_counts():
    # locus is piRNA; its best match is the larger-overlap siRNA region (fp), but the
    # piRNA region it also overlaps (>= 0.5 of the region) is still recalled
    regions = [R(0, 60, "piRNA", name="a"), R(50, 400, "siRNA", name="b")]
    loci = [L(0, 200, "piRNA")]
    m = classification_metrics(loci, regions, CFG)["per_class"]
    assert m["piRNA"]["recall"] == 1.0
    assert m["piRNA"]["tp"] == 0 and m["piRNA"]["fp"] == 1


def test_recall_expressed_counts_any_overlap_over_expressed_only():
    regions = [R(0, 100, "piRNA", name="a"), R(110, 210, "piRNA", name="b"),
               R(220, 320, "piRNA", name="c")]
    expressed = {id(regions[0]), id(regions[1])}
    pc = _pc([L(0, 210, "piRNA")], regions, expressed=expressed)
    assert pc["recall"] == pytest.approx(2 / 3)
    assert pc["recall_expressed"] == pytest.approx(1.0)
