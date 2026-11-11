"""Per-cluster figure data for the HTML report."""
from collections import namedtuple

import bamutil
import class_signals
import report_signals as rs

Aln = namedtuple("Aln", "chrom pos seq nh strand")
Cluster = namedtuple("Cluster", "name chrom start end")


def test_sizes_weighted_by_nh():
    alns = [Aln("c1", 0, "A" * 21, 1, "+"), Aln("c1", 0, "A" * 22, 2, "-")]
    d = rs.cluster_signals(alns, 0, 100)
    assert d["sizes"]["lo"] == 18
    assert d["sizes"]["plus"][21 - 18] == 1.0
    assert d["sizes"]["minus"][22 - 18] == 0.5


def test_nt_position1_is_U_for_T_start():
    alns = [Aln("c1", 0, "T" + "A" * 20, 1, "+") for _ in range(3)]
    d = rs.cluster_signals(alns, 0, 100)
    assert len(d["nt"]) == 26
    assert d["nt"][0]["U"] == 1.0
    assert d["nt"][25] == {"A": 0, "C": 0, "G": 0, "U": 0}


def test_pingpong_peak_at_10():
    alns = [Aln("c1", 100, "A" * 20, 1, "+"), Aln("c1", 90, "A" * 20, 1, "-")]
    d = rs.cluster_signals(alns, 0, 300)
    assert len(d["pingpong"]) == 30
    assert d["pingpong"][9] == max(d["pingpong"]) > 0


def test_coverage_binned():
    alns = [Aln("c1", 10, "A" * 22, 1, "+"), Aln("c1", 5000, "A" * 22, 1, "-")]
    d = rs.cluster_signals(alns, 0, 10000)
    cov = d["coverage"]
    assert len(cov["plus"]) <= 120 and len(cov["minus"]) <= 120
    plus, _ = class_signals.coverage_track(alns, 0, 10000)
    assert round(sum(cov["plus"])) == sum(plus)
    assert cov["start"] == 0 and cov["bin"] >= 1


def test_zero_alignments_zero_filled():
    d = rs.cluster_signals([], 0, 500)
    assert d["n_alns"] == 0
    assert not any(d["sizes"]["plus"]) and not any(d["pingpong"])
    assert len(d["nt"]) == 26 and all(sum(p.values()) == 0 for p in d["nt"])
    assert not any(d["coverage"]["plus"])


def test_missing_chromosome(tmp_path):
    bam = tmp_path / "t.bam"
    bamutil.make_bam(bam, [("c1", 100, False), ("c1", 110, False)])
    out = rs.signals_for(str(bam), [Cluster("a", "chrZ", 0, 500),
                                    Cluster("b", "c1", 0, 500)], xm_as_nh=False)
    assert out["a"]["n_alns"] == 0 and len(out["a"]["nt"]) == 26
    assert out["b"]["n_alns"] == 2
