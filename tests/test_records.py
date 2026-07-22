# tests/test_records.py
from records import ClusterRecord

def test_clusterrecord_defaults_and_length():
    r = ClusterRecord(chrom="chr1", start=100, end=250)
    assert r.length == 150
    assert r.cls == ""
    assert r.review is False
    assert r.dicer_call == "N"
    assert r.size_counts == {}
    assert r.class_scores == {}
