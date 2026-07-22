"""Reading of alignments (SAM/BAM) for the classification engine."""
import os
from collections import namedtuple
import pysam

Aln = namedtuple("Aln", ["chrom", "pos", "strand", "seq", "nh"])


def _keep(rec):
    """Whether an alignment record is a primary mapped read we count."""
    if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
        return False
    return bool(rec.query_sequence)


def _to_aln(rec):
    nh = rec.get_tag("NH") if rec.has_tag("NH") else 1
    strand = "-" if rec.is_reverse else "+"
    return Aln(rec.reference_name, rec.reference_start, strand, rec.query_sequence, nh)


def read_alignments(sam_path):
    """Iterates over primary mapped alignments. Excludes secondary (flag 256)."""
    mode = "rb" if str(sam_path).endswith(".bam") else "r"
    with pysam.AlignmentFile(sam_path, mode) as fh:
        for rec in fh:
            if _keep(rec):
                yield _to_aln(rec)


def read_region(bam_path, chrom, start=None, end=None):
    """Iterates primary mapped alignments in a region of an INDEXED bam (fetch)."""
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        for rec in fh.fetch(chrom, start, end):
            if _keep(rec):
                yield _to_aln(rec)


def ensure_sorted_indexed_bam(sam_path, out_base, threads=1):
    """Return a coordinate-sorted, indexed BAM path for `sam_path`.

    If `sam_path` is already a `.bam` with a `.bai` index, it is returned as-is.
    Otherwise it is sorted (samtools, `-@ threads`) to `<out_base>.sorted.bam`
    and indexed. Sorting is the expensive step and is the one that gets threads.
    """
    sam_path = str(sam_path)
    if sam_path.endswith(".bam") and (
        os.path.exists(sam_path + ".bai") or os.path.exists(sam_path[:-4] + ".bai")
    ):
        return sam_path
    out_bam = f"{out_base}.sorted.bam"
    pysam.sort("-@", str(threads), "-o", out_bam, sam_path)
    pysam.index(out_bam)
    return out_bam


def mapped_chromosomes(bam_path):
    """Reference names with at least one mapped read (uses the BAM index)."""
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        return [s.contig for s in fh.get_index_statistics() if s.mapped > 0]


def total_mapped(bam_path):
    """Total number of mapped alignments (fast, from the BAM index)."""
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        return fh.mapped
