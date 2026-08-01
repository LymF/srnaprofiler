from classplots import Locus, group_by_class, top_n_loci, dominant_dicer


def _l(name, cls, reads, dicer="24"):
    return Locus("chr", 0, 9, name, cls, reads, dicer)


def test_group_by_class_and_restrict():
    loci = [_l("a", "piRNA", 10), _l("b", "siRNA", 5), _l("c", "piRNA", 20)]
    grouped = group_by_class(loci)
    assert set(grouped) == {"piRNA", "siRNA"}
    assert len(grouped["piRNA"]) == 2
    only = group_by_class(loci, classes={"piRNA"})
    assert set(only) == {"piRNA"}


def test_group_by_class_is_case_insensitive():
    # users type "sirna,pirna" but the classifier writes "siRNA"/"piRNA"
    loci = [_l("a", "piRNA", 10), _l("b", "siRNA", 5)]
    grouped = group_by_class(loci, classes={"sirna", "pirna"})
    assert set(grouped) == {"piRNA", "siRNA"}      # canonical casing kept as key
    assert set(group_by_class(loci, classes={"PIRNA"})) == {"piRNA"}


def test_top_n_loci_orders_by_reads():
    loci = [_l("a", "piRNA", 10), _l("c", "piRNA", 20), _l("b", "piRNA", 5)]
    top = top_n_loci(loci, 2)
    assert [l.name for l in top] == ["c", "a"]
    assert top_n_loci(loci, 0) == []


def test_dominant_dicer_ignores_N():
    loci = [_l("a", "piRNA", 1, "24"), _l("b", "piRNA", 1, "24"), _l("c", "piRNA", 1, "N")]
    assert dominant_dicer(loci) == 24
    assert dominant_dicer([_l("a", "piRNA", 1, "N")], default=21) == 21


def test_parse_region():
    from classplots import parse_region
    import pytest
    loc = parse_region("chr1:100-400")
    assert (loc.chrom, loc.start, loc.end) == ("chr1", 100, 400)
    assert loc.name == "chr1:100-400"
    with pytest.raises(ValueError):
        parse_region("chr1-100-400")
