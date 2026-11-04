import gzip
from pathlib import Path

import pytest

from benchmark import fetch
from benchmark.manifest import Library


def make_lib(url: str = "https://example.org/g/genome.fa.gz") -> Library:
    return Library("SRR000001", "Bombyx mori", "test", "arthropod", "ovary", 1, url, ["x"])


def fake_run_factory(calls, files):
    def fake_run(cmd, **kwargs):
        calls.append(cmd)
        out = Path(cmd[cmd.index("-o") + 1])
        out.mkdir(parents=True, exist_ok=True)
        for name in files:
            (out / name).write_bytes(b"data")

        class R:
            returncode = 0

        return R()

    return fake_run


def test_fetch_library_skips_existing(tmp_path, monkeypatch):
    (tmp_path / "SRR000001.fastq.gz").write_bytes(b"x")
    calls = []
    monkeypatch.setattr(fetch.subprocess, "run", fake_run_factory(calls, []))
    assert fetch.fetch_library(make_lib(), tmp_path) == tmp_path / "SRR000001.fastq.gz"
    assert calls == []


def test_fetch_library_calls_iseq_with_accession(tmp_path, monkeypatch):
    calls = []
    monkeypatch.delenv("ISEQ_BIN", raising=False)
    monkeypatch.setattr(fetch.subprocess, "run", fake_run_factory(calls, ["SRR000001.fastq.gz"]))
    out = fetch.fetch_library(make_lib(), tmp_path)
    assert out == tmp_path / "SRR000001.fastq.gz" and out.exists()
    cmd = calls[0]
    assert cmd[0] == "iseq"
    assert cmd[cmd.index("-i") + 1] == "SRR000001"
    assert "-g" in cmd


def test_fetch_library_renames_single_differently_named_file(tmp_path, monkeypatch):
    monkeypatch.setattr(fetch.subprocess, "run", fake_run_factory([], ["SRR000001_xyz.fastq.gz"]))
    out = fetch.fetch_library(make_lib(), tmp_path)
    assert out.name == "SRR000001.fastq.gz" and out.exists()


def test_iseq_bin_env_var(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setenv("ISEQ_BIN", "/opt/iseq")
    monkeypatch.setattr(fetch.subprocess, "run", fake_run_factory(calls, ["SRR000001.fastq.gz"]))
    fetch.fetch_library(make_lib(), tmp_path)
    assert calls[0][0] == "/opt/iseq"


@pytest.mark.parametrize("files", [["SRR000001_1.fastq.gz", "SRR000001_2.fastq.gz"],
                                   ["SRR000001_R1.fastq.gz", "SRR000001_R2.fastq.gz"]])
def test_paired_run_rejected(tmp_path, monkeypatch, files):
    monkeypatch.setattr(fetch.subprocess, "run", fake_run_factory([], files))
    with pytest.raises(ValueError, match="paired"):
        fetch.fetch_library(make_lib(), tmp_path)


def test_fetch_genome_gunzips(tmp_path, monkeypatch):
    def fake_download(url, target):
        with gzip.open(target, "wt") as f:
            f.write(">chr1\nACGT\n")

    monkeypatch.setattr(fetch, "_download", fake_download)
    out = fetch.fetch_genome(make_lib(), tmp_path)
    assert out == tmp_path / "genome.fa"
    assert out.read_text() == ">chr1\nACGT\n"
    assert not (tmp_path / "genome.fa.gz").exists()


def test_fetch_genome_plain_and_skip(tmp_path, monkeypatch):
    seen = []

    def fake_download(url, target):
        seen.append(url)
        Path(target).write_text(">c\nA\n")

    monkeypatch.setattr(fetch, "_download", fake_download)
    lib = make_lib("https://example.org/g/asm.fasta")
    out = fetch.fetch_genome(lib, tmp_path)
    assert out.suffix == ".fa" and out.read_text() == ">c\nA\n"
    fetch.fetch_genome(lib, tmp_path)
    assert len(seen) == 1
