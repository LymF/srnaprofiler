"""Structural RNA is removed from the reads, before anything is clustered.

tRNA, rRNA, snoRNA and snRNA fragments carry 36.5% of the read mass in the
Arabidopsis library and 25.7% in the Drosophila one, and they are not siRNA or
piRNA. Filtering them at the locus stage is too late -- the locus already
exists and its abundance is already wrong. The reference implementation aligns
the reads against an ncRNA database and keeps only what did not align.
"""
from pathlib import Path

import pytest

from mapping import build_ncrna_filter_cmd


def test_keeps_only_reads_that_did_not_align():
    cmd = build_ncrna_filter_cmd("db_prefix", "in.fa", "clean.fa", threads=4)
    assert "--un" in cmd
    assert cmd[cmd.index("--un") + 1] == "clean.fa"


def test_passes_the_database_and_the_reads():
    cmd = build_ncrna_filter_cmd("db_prefix", "in.fa", "clean.fa", threads=4)
    assert "db_prefix" in cmd and "in.fa" in cmd


def test_allows_one_mismatch_like_the_reference():
    cmd = build_ncrna_filter_cmd("db_prefix", "in.fa", "clean.fa", threads=4)
    assert cmd[cmd.index("-v") + 1] == "1"


def test_threads_are_forwarded():
    cmd = build_ncrna_filter_cmd("db_prefix", "in.fa", "clean.fa", threads=8)
    assert cmd[cmd.index("-p") + 1] == "8"


def test_fasta_input_is_declared():
    """bowtie needs telling the input is FASTA, not FASTQ."""
    assert "-f" in build_ncrna_filter_cmd("db", "in.fa", "out.fa", threads=1)


def test_the_alignments_themselves_are_discarded():
    """We want the reads that did NOT match; the SAM of matches is waste."""
    cmd = build_ncrna_filter_cmd("db", "in.fa", "out.fa", threads=1)
    assert "/dev/null" in cmd


def test_fastq_input_is_declared():
    """The library arrives as FASTQ far more often than as FASTA."""
    cmd = build_ncrna_filter_cmd("db", "in.fq", "out.fq", threads=1, fmt="fastq")
    assert "-q" in cmd and "-f" not in cmd


def test_counts_reads_in_fasta(tmp_path):
    """The filter has to report how much it removed: 3% and 90% are not the same."""
    from mapping import count_reads

    fa = tmp_path / "in.fa"
    fa.write_text("".join(f">r{i}\nACGTACGTACGTACGTACGTAC\n" for i in range(10)))
    assert count_reads(str(fa)) == 10


def test_counts_reads_in_fastq(tmp_path):
    from mapping import count_reads

    fq = tmp_path / "in.fq"
    fq.write_text("".join(f"@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n" for i in range(7)))
    assert count_reads(str(fq)) == 7


def test_counts_reads_in_gzipped_fastq(tmp_path):
    import gzip
    from mapping import count_reads

    fq = tmp_path / "in.fq.gz"
    with gzip.open(fq, "wt") as fh:
        fh.write("".join(f"@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n" for i in range(5)))
    assert count_reads(str(fq)) == 5


def _reads(path, n, fasta=True):
    if fasta:
        path.write_text("".join(f">r{i}\nACGTACGTACGTACGTACGTAC\n" for i in range(n)))
    else:
        path.write_text("".join(f"@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n" for i in range(n)))
    return str(path)


def test_filter_returns_the_clean_reads(tmp_path, monkeypatch):
    import mapping

    reads_in = _reads(tmp_path / "in.fa", 10)
    seen = {}

    def fake_run(cmd, env=None):
        seen["cmd"] = cmd
        _reads(Path(cmd[cmd.index("--un") + 1]), 4)

    monkeypatch.setattr(mapping, "run_command", fake_run)
    monkeypatch.setattr(mapping, "check_bowtie_index", lambda p: True)

    out = mapping.filter_ncrna_reads("db", reads_in, tmp_path / "clean.fa", threads=2)

    assert Path(out).exists()
    assert seen["cmd"][0] == "bowtie"
    assert mapping.count_reads(out) == 4


def test_filter_reports_how_much_was_removed(tmp_path, monkeypatch, capsys):
    """A filter that drops 3% and one that drops 90% are not the same run."""
    import mapping

    reads_in = _reads(tmp_path / "in.fa", 10)
    monkeypatch.setattr(mapping, "run_command",
                        lambda cmd, env=None: _reads(Path(cmd[cmd.index("--un") + 1]), 4))
    monkeypatch.setattr(mapping, "check_bowtie_index", lambda p: True)

    mapping.filter_ncrna_reads("db", reads_in, tmp_path / "clean.fa")

    captured = capsys.readouterr()
    logged = captured.out + captured.err
    assert "6 of 10" in logged and "60.0%" in logged


def test_filter_refuses_a_missing_database(tmp_path, monkeypatch):
    """Silently mapping against nothing would leave every read in place."""
    import mapping

    reads_in = _reads(tmp_path / "in.fa", 3)
    monkeypatch.setattr(mapping, "check_bowtie_index", lambda p: False)

    with pytest.raises(RuntimeError, match="ncRNA"):
        mapping.filter_ncrna_reads("db", reads_in, tmp_path / "clean.fa")


def test_filter_keeps_the_input_format(tmp_path, monkeypatch):
    import mapping

    reads_in = _reads(tmp_path / "in.fq", 5, fasta=False)
    seen = {}
    monkeypatch.setattr(mapping, "check_bowtie_index", lambda p: True)
    monkeypatch.setattr(mapping, "run_command",
                        lambda cmd, env=None: (seen.update(cmd=cmd),
                                               _reads(Path(cmd[cmd.index("--un") + 1]), 2,
                                                      fasta=False)))

    mapping.filter_ncrna_reads("db", reads_in, tmp_path / "clean.fq", fmt="fastq")
    assert "-q" in seen["cmd"]
