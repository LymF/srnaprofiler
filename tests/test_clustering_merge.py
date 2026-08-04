"""Bilateral extension, as ShortStack does it (bedtools slop + merge).

Each island grows by `pad` on both sides and overlapping results are merged, so
the effective tolerated gap is 2*pad -- not pad, which is what the old
gap-based rule gave.
"""
from clustering import extend_and_merge


def test_islands_closer_than_two_pads_merge():
    # gap of 150 between them, pad 100 -> each grows 100, so they overlap
    assert extend_and_merge([(1000, 1100), (1250, 1350)], pad=100) == [(900, 1450)]


def test_islands_further_than_two_pads_stay_separate():
    # gap of 250, pad 100 -> extended ranges do not touch
    assert extend_and_merge([(1000, 1100), (1350, 1450)], pad=100) == [(900, 1200), (1250, 1550)]


def test_extension_is_clamped_at_zero_and_chrom_end():
    assert extend_and_merge([(10, 50)], pad=100, chrom_len=200) == [(0, 150)]
    assert extend_and_merge([(150, 190)], pad=100, chrom_len=200) == [(50, 200)]


def test_pad_zero_returns_islands_untouched():
    assert extend_and_merge([(10, 20), (40, 50)], pad=0) == [(10, 20), (40, 50)]


def test_empty_input():
    assert extend_and_merge([], pad=100) == []


def test_touching_islands_merge_into_one():
    assert extend_and_merge([(100, 200), (200, 300)], pad=0) == [(100, 300)]
