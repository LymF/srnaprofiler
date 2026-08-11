"""miRScore's results read as rows, with a missing column failing loud."""
import pytest

from mirna import read_mirscore_rows, parse_mirscore_results

# Header and rows copied verbatim from a real miRScore_results.csv
# (/mnt/d/srnaprofiler-teste/plant/mirna_fresh/mirna/validation/mirscore).
HEADER = ("name,mSeq,mLen,mStart,mStop,msSeq,msLen,msStart,msStop,precSeq,"
          "precLen,result,flags,mismatches,bulges,mReads,msReads,totReads,"
          "numLibraries,precision\n")
ROWS = ("NC_003070.9_1420,CUGAAGUGUUUGGGGGAACUC,21,74.0,94.0,"
        "GUUCCUCUGAGCACUUCAUUG,21.0,24.0,44.0,UUCUUGGUUGGAUGUCUCCUAG,121,"
        "Pass,NA,1.0,0.0,545.0,3.0,550.0,1.0,99.64\n"
        "NC_003074.8_111,CAUGUCGUAAUAGUAGUCACGU,22,25.0,46.0,"
        "GCGAUUUCUAUUAUGACAUAAA,22.0,86.0,107.0,AAGUGAUGAUCAGCUAUCGUAU,133,"
        "Fail,No mature or star reads detected,3.0,0.0,149.0,0.0,164.0,1.0,"
        "90.85\n")


def write_csv(tmp_path, text):
    p = tmp_path / "miRScore_results.csv"
    p.write_text(text)
    return p


def test_rows_are_keyed_by_name_and_keep_every_column(tmp_path):
    rows = read_mirscore_rows(write_csv(tmp_path, HEADER + ROWS))
    assert set(rows) == {"NC_003070.9_1420", "NC_003074.8_111"}
    assert rows["NC_003070.9_1420"]["msSeq"] == "GUUCCUCUGAGCACUUCAUUG"
    assert rows["NC_003074.8_111"]["flags"] == "No mature or star reads detected"
    assert rows["NC_003070.9_1420"]["precision"] == "99.64"


def test_a_missing_required_column_fails_loud(tmp_path):
    """Silently writing NA is how a column rename becomes a wrong table."""
    header = HEADER.replace(",msReads", ",starReads")
    with pytest.raises(RuntimeError, match="msReads"):
        read_mirscore_rows(write_csv(tmp_path, header + ROWS))


def test_passing_names_still_come_out_of_the_same_read(tmp_path):
    got = parse_mirscore_results(write_csv(tmp_path, HEADER + ROWS))
    assert got == {"NC_003070.9_1420"}


from mirna import mirdeep2_novel_rows, parse_mirdeep2_novel

HDR = ["provisional id", "miRDeep2 score",
       "estimated probability that the miRNA candidate is a true positive",
       "rfam alert", "total read count", "mature read count",
       "loop read count", "star read count",
       "significant randfold p-value", "miRBase miRNA",
       "example miRBase miRNA with the same seed", "UCSC browser",
       "NCBI blastn", "consensus mature sequence", "consensus star sequence",
       "consensus precursor sequence", "precursor coordinate"]

NOVEL = "\n".join([
    "some preamble line",
    "novel miRNAs predicted by miRDeep2",
    "\t".join(HDR),
    "\t".join(["NC_003074.8_25650", "5610.4", "", "-", "11002", "3337", "0",
               "7665", "yes", "-", "-", "-", "-",
               "cuuugucuacaauuuuggaaa", "ucccaaauguagacaaagca",
               "cuuugucuacaauuuuggaaaaagugaugacgccauugcucuuucccaaauguagacaaa",
               "NC_003074.8:3366348..3366411:-"]),
    "\t".join(["NC_003070.9_1025", "10.0", "", "rRNA/tRNA", "500", "300", "0",
               "20", "no", "-", "ath-miR166a", "-", "-",
               "uggacugaagggagcucccuu", "agggaugucuuuuugcaccu",
               "uggacugaagggagcucccuuaaaaagggaugucuuuuugcaccu",
               "NC_003070.9:1234..1325:+"]),
    "",
]) + "\n"


def test_raw_novel_rows_keep_the_mirdeep_only_columns(tmp_path):
    p = tmp_path / "result_01_01_2026.csv"
    p.write_text(NOVEL)
    rows = mirdeep2_novel_rows(p)
    assert [r["provisional id"] for r in rows] == ["NC_003074.8_25650",
                                                   "NC_003070.9_1025"]
    assert rows[1]["example miRBase miRNA with the same seed"] == "ath-miR166a"
    assert rows[1]["rfam alert"] == "rRNA/tRNA"
    assert rows[0]["miRDeep2 score"] == "5610.4"


def test_the_empty_probability_column_stays_empty_rather_than_guessed(tmp_path):
    """miRDeep2 leaves it blank without a reference set; that is its answer."""
    p = tmp_path / "result_01_01_2026.csv"
    p.write_text(NOVEL)
    rows = mirdeep2_novel_rows(p)
    key = "estimated probability that the miRNA candidate is a true positive"
    assert rows[0][key] == ""


def test_the_strand_is_read_from_the_precursor_coordinate(tmp_path):
    """It is right there in the coordinate and was being replaced with a dot."""
    p = tmp_path / "result_01_01_2026.csv"
    p.write_text(NOVEL)
    preds = parse_mirdeep2_novel(p)
    assert [x["strand"] for x in preds] == ["-", "+"]


def test_predictions_keep_their_shape(tmp_path):
    p = tmp_path / "result_01_01_2026.csv"
    p.write_text(NOVEL)
    preds = parse_mirdeep2_novel(p)
    assert preds[1] == {"chrom": "NC_003070.9", "strand": "+",
                        "name": "NC_003070.9_1025", "start": 1234, "end": 1325,
                        "mature": "TGGACTGAAGGGAGCTCCCTT"}


def test_a_coordinate_without_a_strand_is_still_read(tmp_path):
    """Older miRDeep2 output omits it; a missing strand is '.', not a crash."""
    p = tmp_path / "result_01_01_2026.csv"
    p.write_text(NOVEL.replace("NC_003070.9:1234..1325:+", "NC_003070.9:1234..1325"))
    assert parse_mirdeep2_novel(p)[1]["strand"] == "."
