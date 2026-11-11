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


def test_missing_chromosome(tmp_path, capsys):
    bam = tmp_path / "t.bam"
    bamutil.make_bam(bam, [("c1", 100, False), ("c1", 110, False)])
    out = rs.signals_for(str(bam), [Cluster("a", "chrZ", 0, 500),
                                    Cluster("b", "c1", 0, 500)], xm_as_nh=False)
    assert out["a"]["n_alns"] == 0 and len(out["a"]["nt"]) == 26
    assert out["b"]["n_alns"] == 2
    assert "WARN" not in capsys.readouterr().err


def test_fetch_error_is_zero_filled_and_warned(tmp_path, capsys, monkeypatch):
    bam = tmp_path / "t.bam"
    bamutil.make_bam(bam, [("c1", 100, False)])

    def boom(*a, **k):
        raise ValueError("bad record")
    monkeypatch.setattr(rs.samio, "read_region_handle", boom)
    out = rs.signals_for(str(bam), [Cluster("a", "c1", 0, 500),
                                    Cluster("b", "c2", 0, 500)], xm_as_nh=False)
    assert out["a"]["n_alns"] == 0 and len(out["a"]["nt"]) == 26
    err = capsys.readouterr().err
    assert "WARN" in err and "2 cluster(s)" in err and "ValueError: bad record" in err


def test_nt_minus_strand_is_reverse_complemented():
    alns = [Aln("c1", 0, "T" * 1 + "A" * 21, 1, "-")]
    d = rs.cluster_signals(alns, 0, 100)
    # revcomp of T+A*21 = T*21+A -> positions 1-2 are U
    assert d["nt"][0]["U"] == 1.0 and d["nt"][1]["U"] == 1.0
    assert d["nt"][21]["A"] == 1.0


def test_coverage_binning_partial_last_bin():
    alns = [Aln("c1", 0, "A" * 20, 1, "+"), Aln("c1", 1000, "A" * 3, 1, "+")]
    d = rs.cluster_signals(alns, 0, 1003, max_bins=100)
    cov = d["coverage"]
    assert cov["bin"] == 11 and len(cov["plus"]) == 92
    assert sum(cov["plus"]) == 23
    assert cov["plus"][-1] == 2  # last bin spans only 2 bases (1001-1002)


def test_none_coordinate_is_zero_filled_and_warned(tmp_path, capsys):
    bam = tmp_path / "t.bam"
    bamutil.make_bam(bam, [("c1", 100, False)])
    out = rs.signals_for(str(bam), [Cluster("a", "c1", None, 500),
                                    Cluster("b", "c1", 0, 500)], xm_as_nh=False)
    assert out["a"]["n_alns"] == 0 and out["a"].get("failed")
    assert out["b"]["n_alns"] == 1 and not out["b"].get("failed")
    assert "1 cluster(s)" in capsys.readouterr().err


def test_nt_dominant_size_inside_window():
    # many 15-nt fragments must not decide the dominant size
    alns = [Aln("c1", 0, "G" * 15, 1, "+") for _ in range(10)]
    alns += [Aln("c1", 0, "T" + "A" * 20, 1, "+") for _ in range(2)]
    d = rs.cluster_signals(alns, 0, 100)
    assert d["nt"][0]["U"] == 1.0


def test_nt_uses_given_dominant_size():
    alns = [Aln("c1", 0, "T" + "A" * 20, 1, "+") for _ in range(5)]
    alns += [Aln("c1", 0, "C" + "A" * 25, 1, "+") for _ in range(2)]
    d = rs.cluster_signals(alns, 0, 100, dominant=26)
    assert d["nt"][0]["C"] == 1.0
