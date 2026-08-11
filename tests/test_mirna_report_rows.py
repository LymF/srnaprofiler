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
