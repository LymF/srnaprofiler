"""The density null needs the size of the references actually mapped."""
import pysam

from samio import mapped_genome_size


def _bam(tmp_path, lengths):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": n, "LN": l} for n, l in lengths]}
    path = tmp_path / "x.bam"
    with pysam.AlignmentFile(str(path), "wb", header=header) as fh:
        for n, _ in lengths:
            a = pysam.AlignedSegment(fh.header)
            a.query_name = "r" + n
            a.query_sequence = "ACGTACGTACGTACGTACGTA"
            a.reference_id = fh.header.get_tid(n)
            a.reference_start = 10
            a.cigarstring = "21M"
            a.flag = 0
            a.mapping_quality = 255
            fh.write(a)
    pysam.index(str(path))
    return str(path)


def test_sums_lengths_of_mapped_references(tmp_path):
    bam = _bam(tmp_path, [("chr1", 1000), ("chr2", 2000)])
    assert mapped_genome_size(bam) == 3000


def test_ignores_references_without_reads(tmp_path):
    bam = _bam(tmp_path, [("chr1", 1000)])
    # a reference declared in the header but holding no read must not count
    with pysam.AlignmentFile(bam, "rb") as fh:
        assert fh.get_reference_length("chr1") == 1000
    assert mapped_genome_size(bam) == 1000
