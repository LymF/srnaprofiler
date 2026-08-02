"""prepare_genome must emit only ACGTN.

Reference genomes carry IUPAC ambiguity codes (TAIR10 has 469 of W/Y/M/K/R/S/D)
and miRDeep2 refuses the file outright: "Sequences are allowed only to comprise
characters [ACGTNacgtn]".
"""
from pathlib import Path

from qc import prepare_genome


def _read(path):
    return Path(path).read_text()


def test_iupac_codes_become_n(tmp_path):
    src = tmp_path / "g.fa"
    src.write_text(">chr1 with description\nacgtWYMKRSD\n")
    out = _read(prepare_genome(str(src), str(tmp_path), "p"))
    assert "ACGTNNNNNNN" in out


def test_uppercases_and_keeps_acgtn(tmp_path):
    src = tmp_path / "g.fa"
    src.write_text(">chr1\nacgtn\n")
    assert "ACGTN" in _read(prepare_genome(str(src), str(tmp_path), "p"))


def test_headers_are_left_alone(tmp_path):
    # standardize_fasta_headers already trimmed them; do not mangle twice
    src = tmp_path / "g.fa"
    src.write_text(">NC_003070.9\nACGT\n")
    assert ">NC_003070.9" in _read(prepare_genome(str(src), str(tmp_path), "p"))
