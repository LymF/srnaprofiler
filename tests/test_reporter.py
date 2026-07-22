import json
from records import ClusterRecord
from reporter import write_tsv, write_gff3, write_run_json

def _rec():
    r = ClusterRecord(chrom="chr1", start=100, end=200, reads=50, dicer_call="21")
    r.cls, r.confidence, r.review = "siRNA", 0.8, False
    return r

def test_write_tsv(tmp_path):
    p = tmp_path / "out.tsv"
    write_tsv([_rec()], str(p))
    lines = p.read_text().splitlines()
    assert lines[0].startswith("chrom\tstart\tend")
    assert "siRNA" in lines[1]

def test_write_gff3(tmp_path):
    p = tmp_path / "out.gff3"
    write_gff3([_rec()], str(p))
    txt = p.read_text()
    assert txt.startswith("##gff-version 3")
    assert "chr1\tsrnaprofiler\tsiRNA\t101\t200" in txt  # GFF 1-based

def test_write_run_json(tmp_path):
    p = tmp_path / "run.json"
    write_run_json({"mincov": 5}, str(p))
    assert json.loads(p.read_text())["mincov"] == 5
