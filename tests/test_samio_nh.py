"""NH (multi-mapper multiplicity) from bowtie1's XM tag.

bowtie1 emits no NH tag but puts the alignment count in XM (capped at -k). We
read XM as NH only when the @PG header shows bowtie1, because bowtie2's XM means
mismatch count.
"""
from samio import read_alignments, prefers_xm_as_nh


def _write_sam(path, pg_line, xm=3):
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:1000",
        pg_line,
        f"r1\t0\tchr1\t101\t255\t21M\t*\t0\t0\tTAGCTAGCTAGCTAGCTAGCT\t*\tXM:i:{xm}",
    ]
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_prefers_xm_for_bowtie1(tmp_path):
    sam = _write_sam(tmp_path / "b1.sam", "@PG\tID:bowtie\tPN:bowtie\tVN:1.3.1\tCL:bowtie -k 5 idx reads")
    assert prefers_xm_as_nh(sam) is True


def test_does_not_prefer_xm_for_bowtie2(tmp_path):
    sam = _write_sam(tmp_path / "b2.sam", "@PG\tID:bowtie2\tPN:bowtie2\tVN:2.5\tCL:bowtie2 -x idx")
    assert prefers_xm_as_nh(sam) is False


def test_nh_read_from_xm_when_enabled(tmp_path):
    sam = _write_sam(tmp_path / "b1.sam", "@PG\tID:bowtie\tPN:bowtie", xm=3)
    (a,) = list(read_alignments(sam, xm_as_nh=True))
    assert a.nh == 3            # XM read as NH
    (b,) = list(read_alignments(sam, xm_as_nh=False))
    assert b.nh == 1            # XM ignored (default) -> unique
