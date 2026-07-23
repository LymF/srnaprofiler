from samio import Aln
from class_signals import pingpong_offsets, phasing_series


def test_pingpong_peak_at_10():
    # + read 5' at 100 (len 20 -> spans 100..119); - read whose 5' end is at 109
    plus = Aln("chr", 100, "+", "A" * 20, 1)
    minus = Aln("chr", 90, "-", "A" * 20, 1)   # 5' end = 90 + 20 - 1 = 109
    counts = pingpong_offsets([plus, minus], hi=30)
    assert counts[10] == 1.0


def test_phasing_series_peaks_at_period():
    # three + reads with 5' ends 21 nt apart -> phasing period 21
    reads = [Aln("chr", p, "+", "A" * 21, 1) for p in (0, 21, 42)]
    counts = phasing_series(reads, "+", maxdist=100)
    assert counts[21] == 2.0
    assert counts[42] == 1.0


def test_phasing_series_strand_filter():
    reads = [Aln("chr", 0, "-", "A" * 21, 1), Aln("chr", 21, "-", "A" * 21, 1)]
    plus = phasing_series(reads, "+", maxdist=100)
    assert sum(plus.values()) == 0.0     # no + reads
