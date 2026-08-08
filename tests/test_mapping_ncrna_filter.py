"""Structural RNA is removed from the reads, before anything is clustered.

tRNA, rRNA, snoRNA and snRNA fragments carry 36.5% of the read mass in the
Arabidopsis library and 25.7% in the Drosophila one, and they are not siRNA or
piRNA. Filtering them at the locus stage is too late -- the locus already
exists and its abundance is already wrong. The reference implementation aligns
the reads against an ncRNA database and keeps only what did not align.
"""
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
