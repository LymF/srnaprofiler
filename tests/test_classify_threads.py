"""Parallel per-chromosome path: threads must not change the result.

Builds a 2-chromosome SAM, sorts+indexes it, and checks that
`extract_by_chrom` produces identical records with 1 thread (serial) and
2 threads (multiprocessing pool), and that the BAM helpers see both chroms.
"""
from samio import ensure_sorted_indexed_bam, mapped_chromosomes, total_mapped
from features import extract_by_chrom

CFG = {"dicermin": 18, "dicermax": 35, "strand_cutoff": 0.8}


def _sam(path):
    header = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:1000", "@SQ\tSN:chr2\tLN:1000"]
    rows = []

    def rec(name, chrom, pos, seq):
        return "\t".join([name, "0", chrom, str(pos), "255", f"{len(seq)}M",
                          "*", "0", "0", seq, "*", "NH:i:1"])

    # chr1: three 21-nt reads (siRNA-ish window); chr2: three 24-nt reads
    rows += [rec("a1", "chr1", 101, "TAGCTAGCTAGCTAGCTAGCT"),
             rec("a2", "chr1", 103, "GAGCTAGCTAGCTAGCTAGCT"),
             rec("a3", "chr1", 105, "TAGCTAGCTAGCTAGCTAGCA")]
    rows += [rec("b1", "chr2", 201, "TAGCTAGCTAGCTAGCTAGCTAGC"),
             rec("b2", "chr2", 203, "GAGCTAGCTAGCTAGCTAGCTAGC"),
             rec("b3", "chr2", 205, "TAGCTAGCTAGCTAGCTAGCTAGA")]
    path.write_text("\n".join(header + rows) + "\n")
    return str(path)


def _summary(records):
    return sorted((r.chrom, r.start, r.end, r.reads, r.dicer_call, r.frac_top)
                  for r in records)


def test_bam_helpers_and_thread_equivalence(tmp_path):
    sam = _sam(tmp_path / "twochrom.sam")
    bam = ensure_sorted_indexed_bam(sam, str(tmp_path / "tc"), threads=2)

    assert bam.endswith(".sorted.bam")
    assert sorted(mapped_chromosomes(bam)) == ["chr1", "chr2"]
    assert total_mapped(bam) == 6

    chroms = mapped_chromosomes(bam)
    serial = extract_by_chrom(bam, chroms, 1_000_000, CFG, mincov=1, pad=5, threads=1)
    parallel = extract_by_chrom(bam, chroms, 1_000_000, CFG, mincov=1, pad=5, threads=2)

    assert len(serial) == 2                       # one island per chromosome
    assert _summary(serial) == _summary(parallel)  # threads must not change results
