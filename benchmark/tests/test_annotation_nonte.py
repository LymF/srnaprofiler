"""Simple repeats, low-complexity regions and structural-RNA repeats are not TE evidence."""
from annotation import load_layers

GFF = ("##gff-version 3\n"
       "chr1\tRepeatMasker\trepeat_region\t1\t100\t.\t+\t.\tName=(TTTTA)n;Class=Simple_repeat\n"
       "chr1\tRepeatMasker\trepeat_region\t201\t300\t.\t+\t.\tName=A-rich;Class=Low_complexity\n"
       "chr1\tRepeatMasker\trepeat_region\t401\t500\t.\t+\t.\tName=5S;Class=rRNA\n"
       "chr1\tRepeatMasker\trepeat_region\t601\t700\t.\t+\t.\tName=ATCOPIA1;Class=LTR/Copia\n")

OUT = ("   SW   perc perc perc  query  position in query  matching  repeat  position in repeat\n"
       "score   div. del. ins.  sequence begin end (left) repeat class/family begin end (left) ID\n\n"
       "  300  1.0 0.0 0.0  chr1  1  100 (900) + (TA)n Simple_repeat 1 100 (0) 1\n"
       "  300  1.0 0.0 0.0  chr1  201  300 (700) + ATCOPIA1 LTR/Copia 1 100 (0) 2\n")


def _te(tmp_path, text, name):
    p = tmp_path / name
    p.write_text(text)
    return load_layers([f"te:{p}"], {"chr1"})[0]


def test_gff_non_te_repeats_left_out(tmp_path):
    L = _te(tmp_path, GFF, "r.gff3")
    assert L.n_features == 1


def test_repeatmasker_out_simple_repeat_left_out(tmp_path):
    L = _te(tmp_path, OUT, "r.out")
    assert L.n_features == 1


def test_bed_class_slash_family_simple_repeat_left_out(tmp_path):
    bed = "chr1\t0\t100\tSimple_repeat/(TA)n\t0\t+\nchr1\t200\t300\tLTR/Copia\t0\t+\n"
    L = _te(tmp_path, bed, "r.bed")
    assert L.n_features == 1


def test_bed_simple_repeat_by_name_only(tmp_path):
    bed = ("chr1\t0\t100\t(TTTTA)n\t0\t+\nchr1\t100\t200\tA-rich\t0\t+\n"
           "chr1\t200\t300\tRNA\t0\t+\nchr1\t400\t500\tATCOPIA1\t0\t+\n")
    L = _te(tmp_path, bed, "n.bed")
    assert L.n_features == 1


def test_generic_rna_class_left_out(tmp_path):
    gff = ("##gff-version 3\n"
           "chr1\tRepeatMasker\trepeat_region\t1\t100\t.\t+\t.\tName=x;Class=RNA\n")
    assert _te(tmp_path, gff, "r.gff3").n_features == 0
