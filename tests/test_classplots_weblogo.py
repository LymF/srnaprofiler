import shutil
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
