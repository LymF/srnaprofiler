"""Regressions found by the first real benchmark run (SRR30544370)."""
import subprocess
from pathlib import Path

import mirna
import qc


def test_filter_fasta_min_length_drops_short_reads(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(
        ">a\n" + "A" * 15 + "\n"
        ">b\n" + "C" * 17 + "\n"
        ">c\n" + "G" * 18 + "\n"
        ">d\n" + "T" * 22 + "\n"
    )
    out = tmp_path / "out.fa"
    kept, dropped = mirna.filter_fasta_min_length(src, out, 18)
    assert (kept, dropped) == (2, 2)
    assert out.read_text() == ">c\n" + "G" * 18 + "\n>d\n" + "T" * 22 + "\n"


def test_filter_fasta_min_length_joins_wrapped_sequences(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">w\n" + "A" * 10 + "\n" + "C" * 10 + "\n>s\nACGT\n")
    out = tmp_path / "out.fa"
    kept, dropped = mirna.filter_fasta_min_length(src, out, 18)
    assert (kept, dropped) == (1, 1)
    assert out.read_text() == ">w\n" + "A" * 10 + "C" * 10 + "\n"


def test_fastqc_uses_outdir_flag(monkeypatch, tmp_path):
    seen = []
    monkeypatch.setattr(qc, "which", lambda name: "/bin/fastqc")
    monkeypatch.setattr(qc, "run_command", lambda cmd, *a, **k: seen.append(cmd))
    qc.fastqc("reads.fq", tmp_path, 2)
    cmd = seen[0]
    assert "--outdir" in cmd and cmd[cmd.index("--outdir") + 1] == str(tmp_path)
    assert "-o" not in cmd


# --- adapter detection, compound read extensions, quantifier without references ---

import gzip

import pytest


def _write_fastq(path, seqs):
    text = "".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(seqs))
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        Path(path).write_text(text)


BIO = "ACGTTGCAGGCTAACCGTTAGC"


def test_detect_adapter_truseq(tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "AGATCGGAAGAGCACAC"] * 10)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS) == "AGATCGGAAGAGC"


def test_detect_adapter_illumina_srna(tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "TGGAATTCTCGGGTGCC"] * 10)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS) == qc.FASTQ_ADAPTERS["illumina_srna"]


def test_detect_adapter_none_when_absent(tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "CCCCCCCCCCCC"] * 10)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS) is None


def test_detect_adapter_plain_file(tmp_path):
    fq = tmp_path / "t.fastq"
    _write_fastq(fq, [BIO + "AGATCGGAAGAGCAC"] * 10)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS) == "AGATCGGAAGAGC"


def test_detect_adapter_respects_read_cap(tmp_path):
    fq = tmp_path / "t.fastq"
    _write_fastq(fq, [BIO + "CCCCCCCCCCCC"] * 20 + [BIO + "AGATCGGAAGAGCAC"] * 20)
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS, n_reads=20) is None
    assert qc.detect_adapter(str(fq), qc.FASTQ_ADAPTERS, n_reads=40) == "AGATCGGAAGAGC"


def _run_fastp(monkeypatch, tmp_path, reads, adapter, in2=None):
    seen = []
    monkeypatch.setattr(qc, "which", lambda name: "/bin/fastp")
    monkeypatch.setattr(qc, "run_command", lambda cmd, *a, **k: seen.append(cmd))
    out2 = tmp_path / "o2.fq" if in2 else None
    qc.trim_with_fastp(str(reads), tmp_path / "o1.fq", in2, out2, 18, 20,
                       adapter, 1, 33, "lib")
    return seen[0]


def test_fastp_detects_truseq_when_no_adapter(monkeypatch, tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "AGATCGGAAGAGCAC"] * 10)
    cmd = _run_fastp(monkeypatch, tmp_path, fq, None)
    assert cmd[cmd.index("--adapter_sequence") + 1] == "AGATCGGAAGAGC"


def test_fastp_single_end_without_adapter_passes_none(monkeypatch, tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "CCCCCCCCCCCC"] * 10)
    cmd = _run_fastp(monkeypatch, tmp_path, fq, None)
    assert "--adapter_sequence" not in cmd


def test_fastp_paired_end_keeps_detect_flag(monkeypatch, tmp_path):
    fq = tmp_path / "t.fastq.gz"
    _write_fastq(fq, [BIO + "CCCCCCCCCCCC"] * 10)
    cmd = _run_fastp(monkeypatch, tmp_path, fq, None, in2=str(fq))
    assert "--detect_adapter_for_pe" in cmd and "--adapter_sequence" not in cmd


def test_fastp_explicit_adapter_skips_detection(monkeypatch, tmp_path):
    def boom(*a, **k):
        raise AssertionError("detect_adapter must not be called")
    monkeypatch.setattr(qc, "detect_adapter", boom)
    cmd = _run_fastp(monkeypatch, tmp_path, tmp_path / "x.fq", "ACGTACGT")
    assert cmd[cmd.index("--adapter_sequence") + 1] == "ACGTACGT"


@pytest.mark.parametrize("name,expected", [
    ("x.input_reads.fastq.gz", "x_input_reads.fastq.gz"),
    ("x.ncrna_filtered.fastq", "x_ncrna_filtered.fastq"),
    ("reads.fq.gz", "reads.fq.gz"),
    ("reads", "reads"),
    ("a.b.fa.bz2", "a_b.fa.bz2"),
])
def test_staged_reads_name_keeps_compound_extension(name, expected):
    assert mirna.staged_reads_name(name) == expected


def test_quantifier_skipped_without_references(monkeypatch, tmp_path):
    def boom(*a, **k):
        raise AssertionError("quantifier.pl must not run")
    monkeypatch.setattr(mirna.subprocess, "run", boom)
    assert mirna.run_quantifier(tmp_path / "r.fa", tmp_path, "Drosophila melanogaster") is True
