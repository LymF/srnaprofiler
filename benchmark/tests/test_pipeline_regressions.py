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
