"""read_region_handle: fetch a region from an already-open BAM handle.

Reusing one open handle across many regions avoids the per-region open/close
overhead that dominates when a run has tens of thousands of loci.
"""
import pysam
from samio import read_region_handle


def _make_bam(path):
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 500}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, pos in enumerate((10, 20, 300)):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.query_sequence = "A" * 20
            a.flag = 0
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 255
            a.cigar = [(0, 20)]
            bam.write(a)
    pysam.sort("-o", str(path), str(path))
    pysam.index(str(path))


def test_read_region_handle_fetches_only_the_region(tmp_path):
    bam = tmp_path / "reads.bam"
    _make_bam(bam)
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        first = [a.pos for a in read_region_handle(fh, "chr1", 0, 100)]
        second = [a.pos for a in read_region_handle(fh, "chr1", 200, 400)]
    assert sorted(first) == [10, 20]      # reads in the first window only
    assert second == [300]                # the same handle is reusable across regions
