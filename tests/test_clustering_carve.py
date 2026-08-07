"""Known miRNA loci are carved out of the clusters that would swallow them.

A miRNA hairpin is ~90 bp and a cluster is ~600-1000 bp, so the hairpin ends up
as 14% of a locus whose class is decided by everything around it. Carving makes
it a locus of its own and leaves the flanks as theirs.
"""
from clustering import carve_out


def test_a_mirna_inside_a_locus_splits_it_in_three():
    assert carve_out([(1000, 2000)], [(1400, 1500)]) == [
        (1000, 1400), (1400, 1500), (1500, 2000)]


def test_a_mirna_at_the_start_leaves_no_empty_piece():
    assert carve_out([(1000, 2000)], [(1000, 1100)]) == [(1000, 1100), (1100, 2000)]


def test_a_mirna_at_the_end_leaves_no_empty_piece():
    assert carve_out([(1000, 2000)], [(1900, 2000)]) == [(1000, 1900), (1900, 2000)]


def test_a_mirna_is_clipped_to_the_locus_it_overlaps():
    # the annotation runs past the locus; only the overlap becomes a locus
    assert carve_out([(1000, 2000)], [(1900, 2500)]) == [(1000, 1900), (1900, 2000)]


def test_a_mirna_covering_the_whole_locus_leaves_it_alone():
    assert carve_out([(1000, 2000)], [(900, 2100)]) == [(1000, 2000)]


def test_two_mirnas_in_one_locus():
    assert carve_out([(1000, 2000)], [(1200, 1300), (1600, 1700)]) == [
        (1000, 1200), (1200, 1300), (1300, 1600), (1600, 1700), (1700, 2000)]


def test_a_mirna_outside_every_locus_changes_nothing():
    assert carve_out([(1000, 2000)], [(5000, 5100)]) == [(1000, 2000)]


def test_slivers_below_the_minimum_are_dropped():
    # a 5-bp flank is not a locus; the miRNA piece itself is always kept
    assert carve_out([(1000, 2000)], [(1005, 1100)], min_piece=20) == [
        (1005, 1100), (1100, 2000)]


def test_output_stays_sorted_and_disjoint():
    out = carve_out([(1000, 2000), (3000, 4000)], [(1500, 1600), (3800, 3900)])
    assert out == sorted(out)
    assert all(a[1] <= b[0] for a, b in zip(out, out[1:]))


def test_no_mirna_intervals_is_a_no_op():
    assert carve_out([(1000, 2000)], []) == [(1000, 2000)]
