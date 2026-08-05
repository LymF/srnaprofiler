"""Bilateral extension: each island grows by `pad` on both sides before
overlaps are merged, so the tolerated gap is 2*pad."""
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


# Merging also requires a compatible profile: distance alone would fuse a 24-nt
# neighbour with a 21-nt one, different pathways.
from clustering import merge_compatible


def _p(size, frac_top):
    return {"size": size, "frac_top": frac_top}


def test_same_profile_and_close_merges():
    islands = [(1000, 1100), (1250, 1350)]
    profiles = [_p(24, 0.5), _p(24, 0.55)]
    assert merge_compatible(islands, profiles, pad=100) == [(900, 1450)]


def test_different_dominant_size_does_not_merge():
    # they split the gap between them: loci must never overlap
    islands = [(1000, 1100), (1250, 1350)]
    profiles = [_p(24, 0.5), _p(21, 0.5)]
    assert merge_compatible(islands, profiles, pad=100) == [(900, 1175), (1175, 1450)]


def test_opposite_strand_bias_does_not_merge():
    islands = [(1000, 1100), (1250, 1350)]
    profiles = [_p(24, 0.05), _p(24, 0.95)]     # difference 0.90 > 0.5
    assert merge_compatible(islands, profiles, pad=100) == [(900, 1175), (1175, 1450)]


def test_unknown_profile_never_blocks_merging():
    # a locus with no size call must not prevent a merge it knows nothing about
    islands = [(1000, 1100), (1250, 1350)]
    profiles = [_p(None, None), _p(24, 0.5)]
    assert merge_compatible(islands, profiles, pad=100) == [(900, 1450)]


def test_chain_merges_only_through_compatible_neighbours():
    # A(24) B(24) C(21): A+B merge, C stays out even though it is near B
    islands = [(1000, 1100), (1250, 1350), (1500, 1600)]
    profiles = [_p(24, 0.5), _p(24, 0.5), _p(21, 0.5)]
    out = merge_compatible(islands, profiles, pad=100)
    assert out == [(900, 1425), (1425, 1700)]


def test_output_is_always_disjoint():
    # reads are assigned to every locus they touch, so overlapping loci would
    # count the same read more than once
    islands = [(1000, 1100), (1150, 1250), (1300, 1400), (1450, 1550)]
    profiles = [_p(21, 0.5), _p(24, 0.5), _p(21, 0.5), _p(24, 0.5)]
    out = merge_compatible(islands, profiles, pad=300)
    assert all(a[1] <= b[0] for a, b in zip(out, out[1:])), out


def test_merge_compatible_empty_input():
    assert merge_compatible([], [], pad=100) == []
