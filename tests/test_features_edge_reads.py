"""FIX 5+6 (TDD): read->region assignment by SPAN OVERLAP.

Coverage in clustering.cluster_regions is built from the read's full span.
Assignment in features.extract must be consistent: a read that starts BEFORE
`start` but enters the island must be counted (previously it was dropped because
only the leftmost position was used).
"""
from features import extract

CONFIG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def test_read_starting_before_region_start_is_counted():
    # edge_read: pos0=90, 20 nt -> span [90,110), starts before start=100.
    # inside_read: pos0=120, 21 nt -> clearly inside.
    region = ("chr1", 100, 160)
    recs = extract("tests/data/edge_reads.sam", [region], lib_size=1_000_000, config=CONFIG)
    r = recs[0]
    # both reads overlap the region by span; the edge_read must NOT be dropped
    assert r.reads == 2
    assert r.distinct_seqs == 2
