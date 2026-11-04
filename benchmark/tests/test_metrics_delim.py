from benchmark.manifest import Config
from benchmark.metrics_delim import delimitation_metrics
from benchmark.truth import Region

CFG = Config(tiny_bp=[50, 100], low_reads=[20, 50, 100], min_overlap_frac=0.5)


def R(start, end, chrom="c1", name="r", cls="piRNA"):
    return Region(chrom, start, end, cls, "t", name)


def L(start, end, reads=500.0, chrom="c1", cls="piRNA"):
    return {"chrom": chrom, "start": start, "end": end, "cls": cls,
            "confidence": 0.9, "reads": reads}


def test_cluster_split_into_three_reports_three_fragments():
    out = delimitation_metrics([L(0, 100), L(300, 400), L(700, 1000)],
                               [R(0, 1000, name="A")], CFG)
    pc = out["per_cluster"][0]
    assert pc["n_fragments"] == 3
    assert pc["recovered_fraction"] == 0.5
    assert (pc["name"], pc["chrom"], pc["start"], pc["end"]) == ("A", "c1", 0, 1000)
    assert out["summary"] == {"n_clusters": 1, "mean_recovered_fraction": 0.5,
                              "mean_fragments": 3.0}


def test_recovered_fraction_is_union_of_overlaps():
    # loci overlap each other; union of region coverage is [100, 400) -> 0.3
    out = delimitation_metrics([L(100, 300), L(200, 400), L(250, 350)],
                               [R(0, 1000)], CFG)
    pc = out["per_cluster"][0]
    assert pc["recovered_fraction"] == 0.3
    assert pc["n_fragments"] == 3


def test_recovered_fraction_clipped_to_region():
    # long locus starting far before the region and ending far after it
    out = delimitation_metrics([L(0, 5000)], [R(1000, 1200)], CFG)
    assert out["per_cluster"][0]["recovered_fraction"] == 1.0
    # long locus starting far before, ending inside the region
    out = delimitation_metrics([L(0, 1100), L(5000, 6000)], [R(1000, 1200)], CFG)
    pc = out["per_cluster"][0]
    assert pc["recovered_fraction"] == 0.5
    assert pc["n_fragments"] == 1


def test_abutting_intervals_do_not_overlap():
    out = delimitation_metrics([L(0, 100), L(200, 300)], [R(100, 200)], CFG)
    pc = out["per_cluster"][0]
    assert pc["n_fragments"] == 0
    assert pc["recovered_fraction"] == 0.0
    assert out["merged_distinct"] == 0


def test_other_chrom_ignored():
    out = delimitation_metrics([L(0, 100, chrom="c2")], [R(0, 100)], CFG)
    assert out["per_cluster"][0]["n_fragments"] == 0


def test_class_is_ignored_including_ambiguous():
    out = delimitation_metrics([L(0, 50, cls="ambiguous"), L(50, 100, cls="miRNA")],
                               [R(0, 100)], CFG)
    pc = out["per_cluster"][0]
    assert pc["n_fragments"] == 2 and pc["recovered_fraction"] == 1.0


def test_locus_spanning_two_clusters_counts_as_merged():
    loci = [L(50, 250),            # spans A and B -> merged
            L(0, 40),              # only A
            L(260, 290),           # only B
            L(100, 200, chrom="c2")]
    clusters = [R(0, 100, name="A"), R(200, 300, name="B"), R(0, 300, chrom="c2", name="C")]
    out = delimitation_metrics(loci, clusters, CFG)
    assert out["merged_distinct"] == 1
    # same locus overlapping 3 clusters is still one locus
    out = delimitation_metrics([L(0, 300)], [R(0, 100), R(100, 200), R(200, 300)], CFG)
    assert out["merged_distinct"] == 1


def test_tiny_counts_use_config_thresholds():
    # lengths 49, 50, 51, 100, 101
    loci = [L(0, 49), L(0, 50), L(0, 51), L(0, 100), L(0, 101)]
    out = delimitation_metrics(loci, [], CFG)
    assert out["tiny_counts"] == {50: 2, 100: 4}
    cfg = Config(tiny_bp=[10], low_reads=[], min_overlap_frac=0.5)
    out = delimitation_metrics([L(0, 10), L(0, 11), L(0, 9)], [], cfg)
    assert out["tiny_counts"] == {10: 2}
    assert out["low_read_counts"] == {}


def test_low_read_counts_strictly_below_threshold():
    loci = [L(0, 10, reads=r) for r in (19.9, 20.0, 49.0, 50.0, 99.0, 100.0, 5000.0)]
    out = delimitation_metrics(loci, [], CFG)
    assert out["low_read_counts"] == {20: 1, 50: 3, 100: 5}


def test_counts_cover_unmatched_loci():
    # loci far from any truth region still count as tiny / low-read
    out = delimitation_metrics([L(9000, 9010, reads=1.0)], [R(0, 100)], CFG)
    assert out["tiny_counts"] == {50: 1, 100: 1}
    assert out["low_read_counts"] == {20: 1, 50: 1, 100: 1}


def test_cluster_order_preserved_and_means():
    out = delimitation_metrics([L(0, 100), L(1000, 1050)],
                               [R(1000, 1100, name="B"), R(0, 100, name="A")], CFG)
    assert [p["name"] for p in out["per_cluster"]] == ["B", "A"]
    assert out["summary"]["n_clusters"] == 2
    assert out["summary"]["mean_recovered_fraction"] == 0.75
    assert out["summary"]["mean_fragments"] == 1.0


def test_empty_inputs():
    out = delimitation_metrics([], [], CFG)
    assert out["per_cluster"] == [] and out["merged_distinct"] == 0
    assert out["tiny_counts"] == {50: 0, 100: 0}
    assert out["low_read_counts"] == {20: 0, 50: 0, 100: 0}
    assert out["summary"] == {"n_clusters": 0, "mean_recovered_fraction": 0.0,
                              "mean_fragments": 0.0}
    out = delimitation_metrics([], [R(0, 100)], CFG)
    assert out["per_cluster"][0]["recovered_fraction"] == 0.0
    assert out["per_cluster"][0]["n_fragments"] == 0
    out = delimitation_metrics([L(0, 10)], [], CFG)
    assert out["tiny_counts"][50] == 1 and out["per_cluster"] == []
