"""Genome browser Task 1: features keep parent and label."""
from annotation import Feature, parse_annotation


def _w(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_old_six_field_construction():
    f = Feature(1, 2, "+", "exon", "")
    assert f.parent == "" and f.label == "" and f.fid == ""


def test_gff3_ensembl_parent_chain_and_label(tmp_path):
    p = _w(tmp_path, "e.gff3",
           "chr1\te\tgene\t1\t100\t.\t+\t.\tID=gene:G1;Name=abc1\n"
           "chr1\te\tmRNA\t1\t100\t.\t+\t.\tID=transcript:T1;Parent=gene:G1\n"
           "chr1\te\texon\t1\t50\t.\t+\t.\tParent=transcript:T1;Name=ex1\n")
    g, m, x = parse_annotation(p)["chr1"]
    assert g.parent == "" and g.label == "abc1"
    assert m.parent == "gene:G1" and m.label == "transcript:T1"
    assert x.parent == "transcript:T1" and x.label == "ex1"


def test_gff3_ncbi_style(tmp_path):
    p = _w(tmp_path, "n.gff3",
           "c\tRefSeq\tgene\t1\t100\t.\t+\t.\tID=gene-LOC1;gene=LOC1\n"
           "c\tRefSeq\tmRNA\t1\t100\t.\t+\t.\tID=rna-XM_1;Parent=gene-LOC1;gene=LOC1\n"
           "c\tRefSeq\tCDS\t5\t90\t.\t+\t0\tID=cds-XP_1;Parent=rna-XM_1\n")
    g, m, c = parse_annotation(p)["c"]
    assert g.label == "LOC1" and g.parent == ""
    assert m.parent == "gene-LOC1"
    assert c.parent == "rna-XM_1" and c.label == "cds-XP_1"


def test_gff3_multiple_parents_take_first(tmp_path):
    p = _w(tmp_path, "m.gff3",
           "c\ts\texon\t1\t10\t.\t+\t.\tParent=t1,t2\n")
    assert parse_annotation(p)["c"][0].parent == "t1"


def test_gtf_parents_and_label(tmp_path):
    p = _w(tmp_path, "a.gtf",
           'c\ts\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "abc";\n'
           'c\ts\ttranscript\t1\t100\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
           'c\ts\texon\t1\t50\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n')
    g, t, x = parse_annotation(p)["c"]
    assert g.parent == "" and g.label == "abc"
    assert t.parent == "G1" and t.label == "G1"
    assert x.parent == "T1" and x.label == "G1"


def test_bed_label(tmp_path):
    p = _w(tmp_path, "a.bed", "c\t0\t10\tpeak1\n" "c\t20\t30\t.\n" "c\t40\t50\n")
    a, b, c = parse_annotation(p)["c"]
    assert a.label == "peak1" and b.label == "" and c.label == ""


def test_rm_out_label_is_repeat_name(tmp_path):
    p = _w(tmp_path, "r.out",
           "  SW  perc perc perc  query  position in query  matching  repeat  class/family  position in repeat  ID\n"
           "  1000 1.0 0.0 0.0 chr1 11 110 (100) + Gypsy-1_DM-I LTR/Gypsy 1 100 (0) 1\n")
    f = parse_annotation(p)["chr1"][0]
    assert f.label == "Gypsy-1_DM-I" and f.name == "LTR/Gypsy"
