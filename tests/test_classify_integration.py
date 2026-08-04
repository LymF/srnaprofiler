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
