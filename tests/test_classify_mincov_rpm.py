"""--mincov is read in reads per million and converted once per run.

As raw counts it did not transfer between libraries: the 10 calibrated on a
30.8M-alignment run means nothing at another depth.
"""
from classify import mincov_to_depth


def test_converts_rpm_using_library_size():
    # 1 RPM of a 10M-read library is 10 reads
    assert mincov_to_depth(1.0, 10_000_000) == 10
    assert mincov_to_depth(0.5, 10_000_000) == 5


def test_same_rpm_gives_different_depth_per_library():
    assert mincov_to_depth(1.0, 5_000_000) == 5
    assert mincov_to_depth(1.0, 40_000_000) == 40


def test_never_below_one():
    # a tiny library must still require at least one read
    assert mincov_to_depth(0.1, 1_000) == 1
    assert mincov_to_depth(0.0, 10_000_000) == 1


def test_handles_unknown_library_size():
    assert mincov_to_depth(1.0, 0) == 1
