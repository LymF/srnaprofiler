"""Finding the genome next to a bowtie index (--index-prefix, no -g).

The glob covered only *.fa and *.fasta, so a genome downloaded from NCBI --
which ships as *.fna -- was not found and the module aborted with "Could not
find genome file from index directory", even though the file was right there.
"""
import pytest

from mirna import find_genome_beside_index


def _touch(d, name):
    p = d / name
    p.write_text(">chr1\nACGT\n")
    return p


def test_finds_ncbi_fna(tmp_path):
    genome = _touch(tmp_path, "GCF_000001735.4_TAIR10.1_genomic.fna")
    (tmp_path / "tair10.1.ebwt").write_text("")
    assert find_genome_beside_index(str(tmp_path / "tair10")) == genome


@pytest.mark.parametrize("name", ["g.fa", "g.fasta", "g.fna", "g.fas"])
def test_accepts_the_usual_fasta_extensions(tmp_path, name):
    genome = _touch(tmp_path, name)
    assert find_genome_beside_index(str(tmp_path / "idx")) == genome


def test_prefers_the_name_matching_the_index(tmp_path):
    # several FASTAs beside the index: picking one at random is a silent wrong
    # answer, so the one sharing the index prefix wins
    _touch(tmp_path, "other.fa")
    mine = _touch(tmp_path, "tair10.fna")
    assert find_genome_beside_index(str(tmp_path / "tair10")) == mine


def test_returns_none_when_absent(tmp_path):
    (tmp_path / "idx.1.ebwt").write_text("")
    assert find_genome_beside_index(str(tmp_path / "idx")) is None


def test_ignores_gzipped_fasta(tmp_path):
    # miRDeep2 needs it uncompressed; silently handing it a .gz fails later and
    # much less clearly
    _touch(tmp_path, "g.fa.gz")
    assert find_genome_beside_index(str(tmp_path / "idx")) is None
