"""Copying reads must keep the format extension.

Path(".fastq.gz").suffix is ".gz", so the copy was named input_reads.gz and the
format detector then rejected it: gzipped FASTQ, the most common input there is,
could not be used.
"""
import pytest

from mirna import reads_suffix
from pathlib import Path


@pytest.mark.parametrize("name,expected", [
    ("s.fastq.gz", ".fastq.gz"),
    ("s.fq.gz", ".fq.gz"),
    ("s.fasta.gz", ".fasta.gz"),
    ("a.trimmed.R1.fastq.gz", ".fastq.gz"),
    ("s.fastq", ".fastq"),
    ("s.fa", ".fa"),
])
def test_keeps_the_format_extension(name, expected):
    assert reads_suffix(Path(name)) == expected


def test_unknown_double_extension_keeps_only_the_last():
    assert reads_suffix(Path("s.weird.gz")) == ".gz"
