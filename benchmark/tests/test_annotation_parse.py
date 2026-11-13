"""Task 1: annotation parsers and interval index."""
import gzip

from annotation import Feature, IntervalIndex, parse_annotation


def _w(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_gff3_exon_and_gene_types(tmp_path):
    p = _w(tmp_path, "a.gff3",
           "##gff-version 3\n"
           "chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\n"
           "chr1\tsrc\texon\t101\t150\t.\t-\t.\tParent=g1\n")
    d = parse_annotation(p)
    assert list(d) == ["chr1"]
    assert d["chr1"][0] == Feature(100, 200, "+", "gene", "gene", "g1")
    assert d["chr1"][1].strand == "-" and d["chr1"][1].ftype == "exon"
    assert d["chr1"][1].start == 100 and d["chr1"][1].end == 150


def test_gff3_te_name_from_classification_attr(tmp_path):
    def name(attrs):
        p = _w(tmp_path, "t.TEanno.gff3",
               f"c\tEDTA\tLTR_retrotransposon\t1\t10\t.\t+\t.\t{attrs}\n")
        return parse_annotation(p)["c"][0].name

    assert name("Name=TE_1;classification=LTR/Gypsy") == "LTR/Gypsy"
    assert name("Name=TE_1;Class=DNA/MITE") == "DNA/MITE"
    assert name("Name=TE_1;class=LINE") == "LINE"
    assert name("Name=TE_1;family=Gypsy") == "Gypsy"
    assert name("Name=TE_1;Family=Copia") == "Copia"
    assert name("Name=TE_1;Class=A;classification=B") == "B"
    assert name("Name=TE_1") == "LTR_retrotransposon"


def test_gff3_ignores_comments_and_fasta_section(tmp_path):
    p = _w(tmp_path, "a.gff",
           "# comment\n"
           "chr1\ts\texon\t1\t10\t.\t+\t.\t.\n"
           "chr1\ts\texon\t50\t20\t.\t+\t.\t.\n"
           "bad line\n"
           "chr1\ts\texon\tx\t20\t.\t+\t.\t.\n"
           "##FASTA\n"
           "chr1\ts\texon\t1\t10\t.\t+\t.\t.\n")
    d = parse_annotation(p)
    assert len(d["chr1"]) == 1


def test_gtf_attributes(tmp_path):
    p = _w(tmp_path, "a.gtf",
           'chr1\ts\texon\t11\t20\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n')
    d = parse_annotation(p)
    f = d["chr1"][0]
    assert f.ftype == "exon" and f.start == 10 and f.end == 20


def test_bed_three_and_six_columns(tmp_path):
    p = _w(tmp_path, "a.bed",
           "track name=x\n#c\nbrowser pos\n"
           "chr1\t5\t15\n"
           "chr1\t20\t30\tfoo\t0\t-\n"
           "chr1\tbad\t30\n")
    d = parse_annotation(p)
    assert d["chr1"][0] == Feature(5, 15, ".", "feature", "")
    assert d["chr1"][1].strand == "-" and d["chr1"][1].start == 20
    assert len(d["chr1"]) == 2


def test_repeatmasker_out(tmp_path):
    hdr = ("   SW  perc perc perc  query  position in query  matching repeat position in repeat\n"
           "score  div. del. ins.  sequence  begin end (left) repeat class/family begin end (left) ID\n"
           "\n")
    r1 = "  500  1.0 0.0 0.0  chr1  101  200  (900)  +  rep1  LTR/Gypsy  1 100 (0)  1\n"
    r2 = "  500  1.0 0.0 0.0  chr1  301  400  (700)  C  rep2  DNA/hAT  (0) 100 1  2\n"
    p = _w(tmp_path, "g.fa.out", hdr + r1 + r2)
    d = parse_annotation(p)
    assert d["chr1"][0] == Feature(100, 200, "+", "feature", "LTR/Gypsy")
    assert d["chr1"][1].strand == "-" and d["chr1"][1].name == "DNA/hAT"
    assert d["chr1"][1].start == 300


def test_gzip_input(tmp_path):
    text = "chr1\ts\tgene\t101\t200\t.\t+\t.\tID=g1\n"
    plain = parse_annotation(_w(tmp_path, "a.gff3", text))
    gz = tmp_path / "b.gff3.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(text)
    assert parse_annotation(str(gz)) == plain
    bed = tmp_path / "c.bed.gz"
    with gzip.open(bed, "wt") as fh:
        fh.write("chr1\t1\t2\n")
    assert parse_annotation(str(bed))["chr1"][0].start == 1


def test_index_overlaps():
    fs = [Feature(0, 10, "+", "a", ""), Feature(5, 100, "+", "b", ""),
          Feature(200, 210, "+", "c", "")]
    ix = IntervalIndex(fs)
    assert ix.overlaps(8, 12) == fs[:2]
    assert ix.overlaps(100, 150) == []
    assert ix.overlaps(205, 206) == [fs[2]]
    big = [Feature(0, 10_000, "+", "big", "")] + [
        Feature(i * 10, i * 10 + 5, "+", "s", "") for i in range(1, 500)]
    assert big[0] in IntervalIndex(big).overlaps(9000, 9001)
    assert IntervalIndex([]).overlaps(0, 5) == []
