"""The evidence scores must reach the TSV, in a fixed column order."""
from records import ClusterRecord
from reporter import write_tsv, _COLS

SCORE_COLS = ["s_1u", "s_5a", "s_fractop", "s_pingpong", "s_phase",
              "s_density", "score"]


def test_score_columns_present_and_after_confidence():
    for c in SCORE_COLS:
        assert c in _COLS
    assert _COLS.index("score") > _COLS.index("confidence")


def test_scores_are_written(tmp_path):
    r = ClusterRecord(chrom="c", start=0, end=100, name="c:0-100", reads=500,
                      dicer_call="24", s_5a=14.2, score=9.4)
    out = tmp_path / "x.tsv"
    write_tsv([r], str(out))
    header, row = out.read_text().splitlines()[:2]
    cols = header.split("\t")
    vals = row.split("\t")
    assert vals[cols.index("s_5a")] == "14.2"
    assert vals[cols.index("score")] == "9.4"


def test_score_defaults_to_zero():
    r = ClusterRecord(chrom="c", start=0, end=1, name="n")
    assert r.score == 0.0
    assert r.s_phase == 0.0
