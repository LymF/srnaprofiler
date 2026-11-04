from bamutil import make_bam
from benchmark.metrics_expr import expressed_regions
from benchmark.truth import Region


def R(chrom, start, end, name):
    return Region(chrom, start, end, "piRNA", "t", name)


def test_expressed_regions_counts_overlapping_alignments(tmp_path):
    bam = tmp_path / "a.bam"
    reads = [("c1", 100 + i, False) for i in range(10)]        # 10 reads in a
    reads += [("c1", 300 + i, False) for i in range(9)]        # 9 reads in b
    reads += [("c1", 500, True)] * 3 + [("c1", 520, False)] * 7  # 7 primary + 3 secondary in c
    make_bam(bam, reads)
    regions = [R("c1", 90, 200, "a"), R("c1", 290, 400, "b"), R("c1", 495, 600, "c"),
               R("c1", 2000, 2100, "empty"), R("missing_contig", 0, 100, "m")]
    got = expressed_regions(bam, regions, min_reads=10)
    names = {r.name for r in regions if id(r) in got}
    assert names == {"a", "c"}  # c counts secondary alignments too; b has 9 < 10


def test_threshold_is_inclusive_and_configurable(tmp_path):
    bam = tmp_path / "a.bam"
    make_bam(bam, [("c2", 100 + i, False) for i in range(5)])
    r = R("c2", 0, 1000, "x")
    assert expressed_regions(bam, [r], min_reads=5) == {id(r)}
    assert expressed_regions(bam, [r], min_reads=6) == set()
