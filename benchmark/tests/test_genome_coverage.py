"""Task 3: window coverage by strand and size class for the genome-browser view."""
import pysam
import pytest

import bamutil
import report_signals as rs
from report_html import dequantise

CHROMS = (("c1", 30000),)
S21, S24, S16 = "A" * 21, "C" * 24, "G" * 16


def _bam(tmp_path, reads):
    bam = tmp_path / "t.bam"
    bamutil.make_seq_bam(bam, reads, CHROMS)
    return str(bam)


def _tracks(bam, chrom, ws, we):
    with pysam.AlignmentFile(bam, "rb") as fh:
        return rs.window_tracks(fh, chrom, ws, we, xm_as_nh=False)


def test_size_classes_constant():
    assert rs.SIZE_CLASSES == ("18-20", "21", "22", "23", "24", "25-35")


def test_depth_and_five_prime_ends(tmp_path):
    bam = _bam(tmp_path, [
        ("c1", 100, S21, "+", 1),    # 100..121
        ("c1", 100, S21, "+", 1),
        ("c1", 200, S24, "-", 1),    # 200..224, 5' end at 223
    ])
    t = _tracks(bam, "c1", 0, 2000)
    assert (t["s"], t["b"], t["nb"]) == (0, 5, 400)
    assert t["k"] == list(rs.SIZE_CLASSES)
    # bins of 5: 100..120 covered by 2 reads over bins 20..24 (bin 24 has 1 base of 5 -> 0.4)
    # quantised to 0..99 of the series max: error <= max / 198
    d21 = dequantise(t["+"]["21"], 400)
    assert d21[20:24] == pytest.approx([2.0] * 4, abs=0.011) and d21[25] == 0
    assert d21[24] == pytest.approx(0.4, abs=0.011)
    d24 = dequantise(t["-"]["24"], 400)
    assert d24[40:44] == pytest.approx([1.0] * 4, abs=0.006)
    assert d24[44] == pytest.approx(0.8, abs=0.006)
    assert t["+"]["24"] == [0, []] and t["-"]["21"] == [0, []]
    assert dequantise(t["e"]["+"], 400)[20] == pytest.approx(2.0)
    en = dequantise(t["e"]["-"], 400)
    assert en[44] == pytest.approx(1.0) and en[40] == 0   # 223 // 5
    assert t["+"]["21"][0] == pytest.approx(2.0)          # per-series max leads the Pairs


def test_nh_weighting(tmp_path):
    bam = _bam(tmp_path, [("c1", 100, S21, "+", 2), ("c1", 100, S21, "+", 2)])
    t = _tracks(bam, "c1", 0, 2000)
    assert t["+"]["21"][0] == pytest.approx(1.0)
    assert t["e"]["+"][0] == pytest.approx(1.0)


def test_outside_window_and_short_reads_excluded(tmp_path):
    bam = _bam(tmp_path, [
        ("c1", 5000, S21, "+", 1),   # outside
        ("c1", 100, S16, "+", 1),    # 16 nt
        ("c1", 100, "A" * 36, "+", 1),   # 36 nt
    ])
    t = _tracks(bam, "c1", 0, 2000)
    for strand in "+-":
        assert all(p == [0, []] for p in t[strand].values())
        assert t["e"][strand] == [0, []]


def test_size_class_edges(tmp_path):
    bam = _bam(tmp_path, [("c1", 100, "A" * n, "+", 1) for n in (18, 20, 22, 23, 25, 35)])
    t = _tracks(bam, "c1", 0, 2000)
    nonzero = {k for k, p in t["+"].items() if p != [0, []]}
    assert nonzero == {"18-20", "22", "23", "25-35"}


def test_bin_size_and_count(tmp_path):
    bam = _bam(tmp_path, [("c1", 100, S21, "+", 1)])
    t = _tracks(bam, "c1", 0, 2000)
    assert (t["b"], t["nb"]) == (5, 400)
    t = _tracks(bam, "c1", 0, 20000)
    assert (t["b"], t["nb"]) == (50, 400)
    t = _tracks(bam, "c1", 10, 1013)
    assert (t["b"], t["nb"]) == (3, 335)


def test_window_offset_and_clipping(tmp_path):
    bam = _bam(tmp_path, [("c1", 95, S21, "+", 1)])   # 95..116, window starts at 100
    t = _tracks(bam, "c1", 100, 120)
    assert (t["s"], t["b"], t["nb"]) == (100, 1, 20)
    assert dequantise(t["+"]["21"], 20) == pytest.approx([1.0] * 16 + [0.0] * 4, abs=0.006)
    assert t["e"]["+"] == [0, []]    # 5' end at 95 is left of the window


def test_for_threads_equal(tmp_path):
    bam = _bam(tmp_path, [("c1", 100 + 7 * i, S21 if i % 2 else S24, "+-"[i % 2], 1 + i % 3)
                          for i in range(40)])
    wins = {f"w{i}": ("c1", i * 50, i * 50 + 1000) for i in range(120)}
    a = rs.window_tracks_for(bam, wins, xm_as_nh=False, threads=1)
    b = rs.window_tracks_for(bam, wins, xm_as_nh=False, threads=2)
    assert a == b and set(a) == set(wins)
    assert (a["w0"]["b"], a["w0"]["nb"]) == (3, 334)
    assert any(p != [0, []] for p in a["w0"]["+"].values())


def test_unknown_contig_zero_filled_one_warn(tmp_path, capsys):
    bam = _bam(tmp_path, [("c1", 100, S21, "+", 1)])
    out = rs.window_tracks_for(bam, {"a": ("chrZ", 0, 2000), "b": ("c1", 0, 2000),
                                     "c": ("chrY", 0, 2000)}, xm_as_nh=False)
    assert out["a"]["+"]["21"] == [0, []] and out["a"]["nb"] == 400
    assert out["b"]["+"]["21"] != [0, []]
    err = capsys.readouterr().err
    assert err.count("WARN") == 1 and "2 window(s)" in err


def test_fetch_error_zero_filled_one_warn(tmp_path, capsys, monkeypatch):
    bam = _bam(tmp_path, [("c1", 100, S21, "+", 1)])

    def boom(*a, **k):
        raise ValueError("bad record")
    monkeypatch.setattr(rs, "_accumulate", boom)
    out = rs.window_tracks_for(bam, {"a": ("c1", 0, 2000)}, xm_as_nh=False)
    assert out["a"]["e"]["+"] == [0, []]
    err = capsys.readouterr().err
    assert err.count("WARN") == 1 and "bad record" in err
