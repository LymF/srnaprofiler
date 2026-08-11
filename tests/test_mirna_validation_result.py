"""Validation carries what was submitted, not only what passed.

miRScore itself is monkeypatched: it is installed in another conda env, so the
suite cannot call it, and a skipif test would skip in the very environment
where the suite runs. Everything around the call is real.
"""
from pathlib import Path

import pytest

import mirna
from mirna import Validation, validate_with_mirscore

#   0         1         2         3         4
#   0123456789012345678901234567890123456789012345
SEQ = "AAAACCCCGGGGTTTTACGTACGTACGTTTTTGGGGCCCCAAAAGG"


class FakeGenome:
    def fetch(self, chrom, start, end):
        return SEQ[start:end]

    def get_reference_length(self, chrom):
        return len(SEQ)


PREDS = [
    # inside the contig with room for flanks: gets re-cut and submitted
    {"chrom": "chr1", "strand": "+", "name": "p1", "start": 17, "end": 24,
     "mature": SEQ[16:24]},
    # mature absent from the genome there: hairpin_with_flank returns None
    {"chrom": "chr1", "strand": "+", "name": "p2", "start": 17, "end": 24,
     "mature": "GGGGGGGG"},
]

# The real header, and a row in the real shape: float-formatted numbers,
# precision on 0-100, sequences in RNA.
RESULTS = ("name,mSeq,mLen,mStart,mStop,msSeq,msLen,msStart,msStop,precSeq,"
           "precLen,result,flags,mismatches,bulges,mReads,msReads,totReads,"
           "numLibraries,precision\n"
           "p1,ACGUACGU,8,9.0,16.0,UUUUGGGG,8.0,25.0,32.0,ACGU,24,Pass,NA,"
           "0.0,0.0,900.0,80.0,1000.0,1.0,99.64\n")


@pytest.fixture
def fake_mirscore(monkeypatch, tmp_path):
    """Stands in for the tool: writes the results file its caller expects."""
    monkeypatch.setattr(mirna.shutil, "which", lambda name: "/fake/miRScore")
    monkeypatch.setattr(mirna, "pysam", None, raising=False)

    def fake_run(cmd, **kwargs):
        out = Path(kwargs["cwd"]) / cmd[cmd.index("-out") + 1]
        out.mkdir(parents=True, exist_ok=True)
        (out / "miRScore_results.csv").write_text(RESULTS)
        return None

    monkeypatch.setattr(mirna.subprocess, "run", fake_run)
    monkeypatch.setattr(mirna, "open_genome", lambda path: FakeGenome(),
                        raising=False)


def test_submitted_holds_the_predictions_that_survived_the_re_cut(fake_mirscore, tmp_path):
    """p2's mature is not in the genome there, so it never reaches miRScore."""
    v = validate_with_mirscore(PREDS, "genome.fa", tmp_path / "reads.fastq",
                               "plant", tmp_path / "work", flank=8)
    assert isinstance(v, Validation)
    assert v.submitted == {"p1"}


def test_rows_carry_mirscores_verdict(fake_mirscore, tmp_path):
    v = validate_with_mirscore(PREDS, "genome.fa", tmp_path / "reads.fastq",
                               "plant", tmp_path / "work", flank=8)
    assert v.rows["p1"]["result"] == "Pass"
    assert v.rows["p1"]["msSeq"] == "UUUUGGGG"   # raw here; DNA happens in build_records


def test_nothing_re_cuttable_means_no_validation_at_all(fake_mirscore, tmp_path):
    """None means 'did not run', which is not the same as 'rejected everything'."""
    bad = [{"chrom": "chr1", "strand": "+", "name": "x", "start": 17, "end": 24,
            "mature": "GGGGGGGG"}]
    assert validate_with_mirscore(bad, "genome.fa", tmp_path / "reads.fastq",
                                  "plant", tmp_path / "work", flank=8) is None


def test_a_missing_miRScore_binary_is_not_validation(monkeypatch, tmp_path):
    monkeypatch.setattr(mirna.shutil, "which", lambda name: None)
    assert validate_with_mirscore(PREDS, "genome.fa", tmp_path / "reads.fastq",
                                  "plant", tmp_path / "work") is None
