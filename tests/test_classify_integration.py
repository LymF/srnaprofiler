import types, os
from classify import run_classify, parse_mirna_bed, mark_mirna
from records import ClusterRecord

def test_mark_mirna_overlap():
    recs = [ClusterRecord("chr1", 100, 200), ClusterRecord("chr1", 500, 600)]
    mark_mirna(recs, {"chr1": [(100, 190)]})   # must dominate the locus, not just overlap
    assert recs[0].cls == "miRNA" and recs[0].review is False
    assert recs[1].cls == ""

def test_run_classify_end_to_end(tmp_path):
    args = types.SimpleNamespace(
        sam="tests/data/toy.sam", genome=None, prefix="t", outdir=str(tmp_path),
        mincov=1, pad=5, mirna_loci=None, config="rules.yaml", norm=1_000_000)
    run_classify(args)
    assert os.path.exists(tmp_path / "t.clusters.tsv")
    assert os.path.exists(tmp_path / "t.clusters.gff3")
    assert os.path.exists(tmp_path / "t.run.json")


def test_a_carved_mirna_locus_reaches_the_tsv_as_unscored(tmp_path):
    """The whole path, not the pieces: mark_mirna -> skipped scorer -> writer.

    The genomic 5' null was wired through every layer except the final call and
    both validations came back byte-identical, so a piece-wise check is not
    enough here either.
    """
    import csv

    bed = tmp_path / "mirna.bed"
    bed.write_text("chr1\t0\t1000\ttest-mir\n")
    args = types.SimpleNamespace(
        sam="tests/data/toy.sam", genome=None, prefix="m", outdir=str(tmp_path),
        mincov=1, pad=5, mirna_loci=str(bed), config="rules.yaml", norm=1_000_000)
    run_classify(args)

    with open(tmp_path / "m.clusters.tsv") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    mirna = [r for r in rows if r["cls"] == "miRNA"]
    assert mirna, "the BED covers the toy locus, so something must be called miRNA"
    for r in mirna:
        assert r["score"] == "NA"
        assert r["confidence"] == "1.0"
    for r in rows:
        if r["cls"] != "miRNA":
            assert r["score"] != "NA"
