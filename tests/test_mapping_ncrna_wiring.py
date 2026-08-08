"""The ncRNA filter has to be reached by run_mapping, not merely exist.

A previous feature (the genomic 5' null) was passed through every layer except
the final call and stayed switched off while both validations came back
byte-identical. What the aligner receives is the only thing that proves the
filter ran.
"""
from pathlib import Path
from types import SimpleNamespace

import pytest

import mapping


@pytest.fixture
def wired(monkeypatch, tmp_path):
    """run_mapping with everything external stubbed; records what bowtie got."""
    seen = {}

    monkeypatch.setattr(mapping, "validate_arguments", lambda args: None)
    monkeypatch.setattr(mapping, "check_disk_space", lambda *a, **k: True)
    monkeypatch.setattr(mapping, "validate_file_exists", lambda *a, **k: None)
    monkeypatch.setattr(mapping, "get_read_stats", lambda p: {"has_uracil": False})
    monkeypatch.setattr(mapping, "check_bowtie_index", lambda p: True)
    monkeypatch.setattr(mapping, "generate_mapping_report", lambda *a, **k: None)

    def fake_align(**kwargs):
        seen["reads1"] = kwargs["reads1"]
        return {}

    monkeypatch.setattr(mapping, "run_bowtie_alignment", fake_align)
    return seen


def _args(tmp_path, **extra):
    reads = tmp_path / "reads.fa"
    reads.write_text(">r1\nACGTACGTACGTACGTACGTAC\n")
    args = SimpleNamespace(
        outdir=str(tmp_path / "out"), prefix="p", reads=str(reads),
        reads_paired=None, genome=None, index_prefix="idx", fmt="auto",
        mismatches=1, k=50, keep_unmapped=False, sam_out=None, threads=1,
        ncrna_db=None,
    )
    for key, value in extra.items():
        setattr(args, key, value)
    return args


def test_filtered_reads_are_what_gets_aligned(wired, tmp_path, monkeypatch):
    clean = tmp_path / "clean.fa"
    clean.write_text(">r1\nACGTACGTACGTACGTACGTAC\n")
    monkeypatch.setattr(mapping, "filter_ncrna_reads",
                        lambda *a, **k: str(clean))

    mapping.run_mapping(_args(tmp_path, ncrna_db="ncrna_idx"))

    assert wired["reads1"] == str(clean)


def test_the_database_and_format_reach_the_filter(wired, tmp_path, monkeypatch):
    seen = {}

    def fake_filter(db_prefix, reads_in, reads_out, threads=1, fmt="fasta"):
        seen.update(db=db_prefix, reads_in=reads_in, fmt=fmt, threads=threads)
        Path(reads_out).parent.mkdir(parents=True, exist_ok=True)
        Path(reads_out).write_text(">r1\nACGTACGTACGTACGTACGTAC\n")
        return str(reads_out)

    monkeypatch.setattr(mapping, "filter_ncrna_reads", fake_filter)

    args = _args(tmp_path, ncrna_db="ncrna_idx", threads=4)
    mapping.run_mapping(args)

    assert seen["db"] == "ncrna_idx"
    assert seen["reads_in"] == args.reads
    assert seen["fmt"] == "fasta"
    assert seen["threads"] == 4


def test_without_the_flag_nothing_is_filtered(wired, tmp_path, monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("filter ran without --ncrna-db")

    monkeypatch.setattr(mapping, "filter_ncrna_reads", refuse)

    args = _args(tmp_path)
    mapping.run_mapping(args)

    assert wired["reads1"] == args.reads
