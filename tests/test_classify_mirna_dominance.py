"""A miRNA only takes a locus's label if it covers most of it.

With bigger loci a 21-nt miRNA can sit inside a cluster holding much more, and
labelling the whole thing miRNA hides everything else.
"""
from records import ClusterRecord
from classify import mark_mirna


def _rec(start, end):
    return ClusterRecord(chrom="c", start=start, end=end,
                         name=f"c:{start}-{end}", reads=100, cls="hcsiRNA")


def test_mirna_covering_the_locus_claims_it():
    r = _rec(1000, 1100)
    mark_mirna([r], {"c": [(1000, 1090)]})
    assert r.cls == "miRNA"


def test_mirna_inside_a_large_locus_does_not_claim_it():
    r = _rec(1000, 11000)                      # 10 kb locus
    mark_mirna([r], {"c": [(5000, 5021)]})     # 21 nt inside it
    assert r.cls == "hcsiRNA"


def test_partial_overlap_below_threshold_does_not_claim():
    r = _rec(1000, 2000)
    mark_mirna([r], {"c": [(1000, 1300)]})     # 30% of the locus
    assert r.cls == "hcsiRNA"


def test_threshold_is_configurable():
    r = _rec(1000, 2000)
    mark_mirna([r], {"c": [(1000, 1300)]}, min_fraction=0.25)
    assert r.cls == "miRNA"


def test_claimed_locus_keeps_the_old_confidence_behaviour():
    r = _rec(1000, 1100)
    mark_mirna([r], {"c": [(1000, 1090)]})
    assert r.confidence == 1.0
    assert r.review is False
