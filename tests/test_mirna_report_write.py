"""The table has everything; the BED and the FASTAs have exactly the passes."""
from mirna_report import (build_records, write_table, write_bed, write_fasta,
                          validated, MIRNA_COLS)

PREDS = [
    {"chrom": "Chr1", "strand": "+", "name": "p1", "start": 1234, "end": 1325,
     "mature": "TGGACTGAAGGGAGCTCCCTT"},
    {"chrom": "Chr2", "strand": "-", "name": "p2", "start": 88100, "end": 88190,
     "mature": "TTGACAGAAGATAGAGAGCAC"},
]
ROWS = {
    "p1": {"result": "Pass", "msSeq": "GUUCCUCUGAGCACUUCAUUG",
           "precSeq": "UUCUUGGUUGGAUGUCUCCUAG", "flags": "NA"},
    "p2": {"result": "Fail", "msSeq": "GCGAUUUCUAUUAUGACAUAAA",
           "precSeq": "AAGUGAUGAUCAGCUAUCGUAU",
           "flags": "No mature or star reads detected"},
}


def recs():
    return build_records(PREDS, "mirdeep-p2", ROWS, submitted={"p1", "p2"})


def test_the_table_has_a_row_per_prediction_and_the_declared_header(tmp_path):
    p = tmp_path / "at.mirna.tsv"
    write_table(recs(), p)
    lines = p.read_text().splitlines()
    assert lines[0].split("\t") == MIRNA_COLS
    assert len(lines) == 3
    assert lines[1].split("\t")[MIRNA_COLS.index("result")] == "Pass"


def test_a_tab_in_a_field_cannot_break_a_row(tmp_path):
    r = recs()
    r[1].flags = "star\tnot\tdetected"
    p = tmp_path / "at.mirna.tsv"
    write_table(r, p)
    for line in p.read_text().splitlines():
        assert len(line.split("\t")) == len(MIRNA_COLS)


def test_the_bed_carries_only_the_passes_in_half_open_coordinates(tmp_path):
    p = tmp_path / "at.mirna_loci.bed"
    assert write_bed(recs(), p) == 1
    assert p.read_text() == "Chr1\t1233\t1325\tp1\n"


def test_the_three_fastas_and_the_bed_are_the_same_set(tmp_path):
    r = recs()
    names = {rec.name for rec in validated(r)}
    for field in ("mature_seq", "star_seq", "precursor_seq"):
        path = tmp_path / f"at.{field}.fa"
        assert write_fasta(r, path, field) == len(names)
        got = {ln[1:] for ln in path.read_text().splitlines() if ln.startswith(">")}
        assert got == names


def test_the_fasta_holds_the_sequence_of_the_requested_field(tmp_path):
    """And in DNA -- miRScore answered in RNA."""
    p = tmp_path / "at.star.fa"
    write_fasta(recs(), p, "star_seq")
    assert p.read_text() == ">p1\nGTTCCTCTGAGCACTTCATTG\n"
