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
