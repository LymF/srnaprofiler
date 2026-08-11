"""One record per prediction: the verdict, the evidence, and the coordinates."""
import pytest

from mirna_report import build_records, MirnaRecord, MIRNA_COLS

PREDS = [
    {"chrom": "Chr1", "strand": "+", "name": "p1", "start": 1234, "end": 1325,
     "mature": "TGGACTGAAGGGAGCTCCCTT"},
    {"chrom": "Chr2", "strand": "-", "name": "p2", "start": 88100, "end": 88190,
     "mature": "TTGACAGAAGATAGAGAGCAC"},
    {"chrom": "Chr3", "strand": "+", "name": "p3", "start": 551020, "end": 551108,
     "mature": "TAAGGGATCTTACAGGTAACG"},
]

# Values copied from a real miRScore_results.csv: RNA sequences, numeric fields
# as float-formatted strings, precision on a 0-100 scale.
SCORE_ROWS = {
    "p1": {"name": "p1", "mSeq": "CUGAAGUGUUUGGGGGAACUC", "mLen": "21",
           "msSeq": "GUUCCUCUGAGCACUUCAUUG", "msLen": "21.0",
           "precSeq": "UUCUUGGUUGGAUGUCUCCUAGAGUUCCUCUGAGCACUUCAUUG",
           "precLen": "121", "mStart": "74.0", "mStop": "94.0",
           "msStart": "24.0", "msStop": "44.0",
           "result": "Pass", "flags": "NA", "mismatches": "1.0", "bulges": "0.0",
           "mReads": "545.0", "msReads": "3.0", "totReads": "550.0",
           "precision": "99.64"},
    "p2": {"name": "p2", "mSeq": "CAUGUCGUAAUAGUAGUCACGU", "mLen": "22",
           "msSeq": "GCGAUUUCUAUUAUGACAUAAA", "msLen": "22.0",
           "precSeq": "AAGUGAUGAUCAGCUAUCGUAUUUCAUGUCGUAAUAGUAGUCACGU",
           "precLen": "133", "mStart": "25.0", "mStop": "46.0",
           "msStart": "86.0", "msStop": "107.0",
           "result": "Fail", "flags": "No mature or star reads detected",
           "mismatches": "3.0", "bulges": "0.0", "mReads": "149.0",
           "msReads": "0.0", "totReads": "164.0", "precision": "90.85"},
}


def test_one_record_per_prediction():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert [r.name for r in recs] == ["p1", "p2", "p3"]


def test_coordinates_become_zero_based_half_open():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert (recs[0].start, recs[0].end) == (1233, 1325)


def test_the_four_verdicts():
    """p1 passed, p2 failed, p3 was never submitted, and nothing is blank."""
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert [r.result for r in recs] == ["Pass", "Fail", "NotScored"]
    assert recs[2].flags == "precursor could not be re-cut from the genome"

    none_ran = build_records(PREDS, "mirdeep-p2", mirscore_rows=None)
    assert {r.result for r in none_ran} == {"NotValidated"}
    assert all(r.star_seq == "NA" for r in none_ran)


def test_submitted_but_absent_from_the_results_is_not_a_failure():
    """miRScore dropping a locus is not miRScore rejecting it."""
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS,
                         submitted={"p1", "p2", "p3"})
    assert recs[2].result == "NotScored"
    assert recs[2].flags == "absent from miRScore results"


def test_submitted_none_means_all_were_submitted():
    """No `submitted` set given: p3 missing from the rows is still absent
    from the results, not blamed on a re-cut failure that never happened."""
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS)
    assert recs[2].result == "NotScored"
    assert recs[2].flags == "absent from miRScore results"


def test_evidence_columns_come_from_mirscore():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert recs[0].star_len == "21.0"        # miRScore's own float string, kept
    assert recs[0].mature_reads == "545.0"
    assert recs[0].precision == "99.64"      # 0-100, not 0-1
    assert recs[0].mature_len == 21          # ours, from the prediction


def test_mirscore_sequences_are_converted_to_dna():
    """miRScore answers in RNA; the mature comes from the prediction in DNA.

    Mixing them would give a mature.fa in T and a star.fa in U.
    """
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert recs[0].star_seq == "GTTCCTCTGAGCACTTCATTG"
    assert "U" not in recs[0].precursor_seq


def test_mirdeep_only_columns_are_na_for_the_plant_engine():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    assert recs[0].mirdeep_score == "NA"
    assert recs[0].mirbase_hit == "NA"
    assert recs[0].engine == "mirdeep-p2"


def test_mirdeep_only_columns_are_filled_when_the_engine_provides_them():
    rows = {"p1": {"provisional id": "p1", "miRDeep2 score": "96.3",
                   "rfam alert": "-", "significant randfold p-value": "yes",
                   "miRBase miRNA": "-",
                   "example miRBase miRNA with the same seed": "ath-miR166a"}}
    recs = build_records(PREDS, "mirdeep2", SCORE_ROWS, submitted={"p1", "p2"},
                         mirdeep_rows=rows)
    assert recs[0].mirdeep_score == "96.3"
    assert recs[0].mirbase_same_seed == "ath-miR166a"
    assert recs[1].mirdeep_score == "NA"     # absent row, not an error


def test_a_repeated_name_fails_loud():
    """The join is by name; a duplicate would silently multiply rows."""
    dup = PREDS + [dict(PREDS[0])]
    with pytest.raises(RuntimeError, match="p1"):
        build_records(dup, "mirdeep-p2", SCORE_ROWS, submitted={"p1"})


def test_every_column_exists_on_the_record():
    recs = build_records(PREDS, "mirdeep-p2", SCORE_ROWS, submitted={"p1", "p2"})
    for col in MIRNA_COLS:
        assert hasattr(recs[0], col), col
