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


def test_a_mirna_inside_a_big_cluster_is_carved_into_its_own_locus(tmp_path):
    """End to end: the carve must reach the clustering, not just exist.

    Without it the hairpin is ~7% of the locus around it, never reaches the
    dominance threshold, and the cluster keeps a siRNA/piRNA label -- which is
    how the miRNA route broke when loci grew from ~100 bp to ~1 kb.
    """
    import types
    from samio import Aln
    import pysam

    # one 1.2 kb cluster of 21-nt reads, with a 90-bp miRNA hairpin inside it
    sam = tmp_path / "in.sam"
    with open(sam, "w") as fh:
        fh.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:c\tLN:5000\n")
        fh.write("@PG\tID:x\tPN:other\n")
        n = 0
        for pos in range(1000, 2200, 12):
            for _ in range(6):
                n += 1
                fh.write(f"r{n}\t0\tc\t{pos}\t255\t21M\t*\t0\t0\t{'A'*21}\t{'I'*21}\tXM:i:1\n")
        for pos in range(1500, 1590, 3):
            for _ in range(6):
                n += 1
                fh.write(f"m{n}\t0\tc\t{pos}\t255\t21M\t*\t0\t0\t{'T'*21}\t{'I'*21}\tXM:i:1\n")

    bed = tmp_path / "mir.bed"
    bed.write_text("c\t1499\t1589\tMIR-test\n")

    from classify import run_classify
    args = types.SimpleNamespace(
        sam=str(sam), genome=None, prefix="m", outdir=str(tmp_path),
        mincov=0.0, pad=200, mirna_loci=str(bed), config="rules.yaml",
        norm=1_000_000, threads=1)
    run_classify(args)

    import csv
    rows = list(csv.DictReader(open(tmp_path / "m.clusters.tsv"), delimiter="\t"))
    mir = [r for r in rows if r["cls"] == "miRNA"]
    assert len(mir) == 1, [(r["start"], r["end"], r["cls"]) for r in rows]
    assert int(mir[0]["start"]) == 1499 and int(mir[0]["end"]) == 1589
    # the flanks survive as their own loci instead of being thrown away
    assert len(rows) >= 3


def test_a_carved_mirna_locus_is_marked_unscored():
    """It never reaches the scorer, so it has no score -- not a score of zero."""
    from classify import mark_mirna
    from records import ClusterRecord

    r = ClusterRecord(chrom="chr1", start=100, end=200)
    mark_mirna([r], {"chr1": [(100, 200)]})
    assert r.cls == "miRNA"
    assert r.score is None


def test_a_locus_the_mirna_does_not_cover_keeps_its_score():
    from classify import mark_mirna
    from records import ClusterRecord

    r = ClusterRecord(chrom="chr1", start=100, end=1100)
    r.score = 0.7
    mark_mirna([r], {"chr1": [(100, 150)]})
    assert r.cls != "miRNA"
    assert r.score == 0.7
