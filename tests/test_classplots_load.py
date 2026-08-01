import pytest
from classplots import load_clusters, Locus

TSV = "tests/data/mini.clusters.tsv"


def test_load_clusters_parses_and_skips_blank_cls():
    loci = load_clusters(TSV)
    assert len(loci) == 2                       # blank-cls row skipped
    assert isinstance(loci[0], Locus)
    piRNA = [l for l in loci if l.cls == "piRNA"][0]
    assert piRNA.chrom == "chr1" and piRNA.start == 100 and piRNA.end == 400
    assert piRNA.reads == 500.0 and piRNA.dicer_call == "24"


def test_load_clusters_missing_column_raises(tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("chrom\tstart\tend\n1\t0\t9\n")
    with pytest.raises(ValueError):
        load_clusters(str(bad))
