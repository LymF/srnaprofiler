from collections import Counter
from signals import offset_spectrum


def test_offset_spectrum_peak_at_10():
    plus5 = {100: 2.0}
    minus5 = {109: 3.0}          # 100 + 10 - 1 = 109 -> overlap length 10
    counts = offset_spectrum(plus5, minus5, lo=1, hi=20)
    assert counts[10] == 6.0
    assert sum(counts.values()) == 6.0


def test_offset_spectrum_same_strand_distance():
    a = {10: 1.0, 31: 1.0, 52: 1.0}   # phased ladder, period 21
    counts = offset_spectrum(a, a, lo=1, hi=60)
    # d is p+d-1 = q, i.e. d = (q - p) + 1, so a period-21 ladder peaks at d=22/43,
    # not 21/42 (the +1 keeps d consistent with the "overlap length" semantics used
    # by the cross-strand ping-pong case, per the mandated p + d - 1 formula).
    assert counts[22] == 2.0          # (10,31) and (31,52)
    assert counts[43] == 1.0          # (10,52)
    assert counts.get(0, 0) == 0      # self-pairs excluded (lo=1)


def test_offset_spectrum_respects_hi():
    counts = offset_spectrum({0: 1.0}, {29: 1.0}, lo=1, hi=20)
    assert counts.get(30, 0) == 0     # overlap length 30 is beyond hi
