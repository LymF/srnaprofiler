import shutil
import subprocess
import classplots
from classplots import run_weblogo


def test_run_weblogo_absent_binary_returns_false(tmp_path, monkeypatch):
    monkeypatch.setattr(shutil, "which", lambda name: None)
    fasta = tmp_path / "s.fa"
    fasta.write_text(">seq_1\nUACG\n")
    assert run_weblogo(str(fasta), str(tmp_path / "o.pdf"), "sense") is False


def test_run_weblogo_empty_fasta_returns_false(tmp_path):
    fasta = tmp_path / "empty.fa"
    fasta.write_text("")
    assert run_weblogo(str(fasta), str(tmp_path / "o.pdf"), "sense") is False


def test_run_weblogo_passes_format_flag(tmp_path, monkeypatch):
    """The `fmt` arg drives `weblogo -F <fmt>`; default is svg (no Ghostscript)."""
    fasta = tmp_path / "s.fa"
    fasta.write_text(">seq_1\nUACG\n")
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/weblogo")
    captured = {}

    def fake_run(cmd, **kw):
        captured["cmd"] = cmd
        class R:
            returncode = 0
            stderr = ""
        return R()

    monkeypatch.setattr(subprocess, "run", fake_run)

    assert run_weblogo(str(fasta), str(tmp_path / "o.svg"), "sense") is True
    cmd = captured["cmd"]
    assert cmd[cmd.index("-F") + 1] == "svg"      # default format

    assert run_weblogo(str(fasta), str(tmp_path / "o.eps"), "sense", fmt="eps") is True
    cmd = captured["cmd"]
    assert cmd[cmd.index("-F") + 1] == "eps"
