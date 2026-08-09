"""A class that bypasses the scorer must not report a score of zero.

miRNA loci are carved out of `--mirna-loci` and never reach the RuleScorer
(classify.py:194), so `score` kept its dataclass default of 0.0. In the
Arabidopsis run that is 294 curated miRBase loci reporting the same 0.0 as the
loci with genuinely no evidence -- sorting the table by score buries every
validated miRNA at the bottom. Zero is a measurement; "not scored" is not.
"""
from records import ClusterRecord
from reporter import write_gff3, write_tsv


def _mirna():
    r = ClusterRecord(chrom="chr1", start=100, end=200, reads=50, dicer_call="21")
    r.cls, r.confidence, r.review, r.score = "miRNA", 1.0, False, None
    return r


def _sirna():
    r = ClusterRecord(chrom="chr1", start=300, end=400, reads=50, dicer_call="21")
    r.cls, r.confidence, r.review, r.score = "siRNA", 0.8, False, 4.2
    return r


def test_an_unscored_locus_reports_na_in_the_tsv(tmp_path):
    p = tmp_path / "out.tsv"
    write_tsv([_mirna()], str(p))
    header, row = p.read_text().splitlines()
    assert row.split("\t")[header.split("\t").index("score")] == "NA"


def test_a_scored_locus_still_reports_its_number(tmp_path):
    p = tmp_path / "out.tsv"
    write_tsv([_sirna()], str(p))
    header, row = p.read_text().splitlines()
    assert row.split("\t")[header.split("\t").index("score")] == "4.2"


def test_a_real_zero_is_not_turned_into_na(tmp_path):
    """Zero is a measured verdict -- no phase evidence at all -- and must stay."""
    r = _sirna()
    r.score = 0.0
    p = tmp_path / "out.tsv"
    write_tsv([r], str(p))
    header, row = p.read_text().splitlines()
    assert row.split("\t")[header.split("\t").index("score")] == "0.0"


def test_the_gff_score_column_uses_the_missing_value_marker(tmp_path):
    """GFF3 column 6 takes a float or '.', never a string like NA."""
    p = tmp_path / "out.gff3"
    write_gff3([_mirna()], str(p))
    row = [ln for ln in p.read_text().splitlines() if not ln.startswith("#")][0]
    assert row.split("\t")[5] == "."
    assert "score=NA" in row
