"""Weighted library size (sum of 1/NH) for RPM.

`total_mapped` counts alignment lines, which with bowtie1 -k5 overcounts a
multi-mapping read up to 5 times. RPM must use the same 1/NH scale as the
per-locus read counts, so the library size is the sum of the weights.
"""
from samio import weighted_library_size


def _write_sam(path, xm_by_read):
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:1000",
        "@PG\tID:bowtie\tPN:bowtie\tVN:1.3.1\tCL:bowtie -k 5 idx reads",
    ]
    for i, xm in enumerate(xm_by_read):
        lines.append(f"r{i}\t0\tchr1\t{101 + i}\t255\t21M\t*\t0\t0\t"
                     f"TAGCTAGCTAGCTAGCTAGCT\t*\tXM:i:{xm}")
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def test_sums_inverse_nh(tmp_path):
    # four alignment lines of a read reported 4 times, plus one unique read
    sam = _write_sam(tmp_path / "w.sam", [4, 4, 4, 4, 1])
    assert weighted_library_size(sam, xm_as_nh=True) == 2.0


def test_ignores_xm_when_not_bowtie1(tmp_path):
    sam = _write_sam(tmp_path / "w.sam", [4, 4, 4, 4, 1])
    assert weighted_library_size(sam, xm_as_nh=False) == 5.0


def test_skips_secondary_alignments(tmp_path):
    sam = _write_sam(tmp_path / "w.sam", [1])
    with open(sam, "a") as fh:
        fh.write("r9\t256\tchr1\t201\t255\t21M\t*\t0\t0\t"
                 "TAGCTAGCTAGCTAGCTAGCT\t*\tXM:i:1\n")
    assert weighted_library_size(sam, xm_as_nh=True) == 1.0
