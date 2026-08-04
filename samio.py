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


def _nh(rec, xm_as_nh):
    """Read multiplicity (number of reported alignments) for 1/NH weighting.

    Prefers the standard NH tag. bowtie1 emits no NH but encodes the alignment
    count in XM (capped at -k), so when the input is bowtie1 (`xm_as_nh`) we use
    XM. We do NOT read XM otherwise: in bowtie2/others XM means mismatch count.
    Defaults to 1 (unique).
    """
    if rec.has_tag("NH"):
        return rec.get_tag("NH")
    if xm_as_nh and rec.has_tag("XM"):
        return max(rec.get_tag("XM"), 1)
    return 1


def _to_aln(rec, xm_as_nh=False):
    strand = "-" if rec.is_reverse else "+"
    return Aln(rec.reference_name, rec.reference_start, strand,
               rec.query_sequence, _nh(rec, xm_as_nh))


def prefers_xm_as_nh(bam_path):
    """True if the @PG header shows bowtie1 (so XM should be read as NH)."""
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        for pg in fh.header.get("PG", []):
            prog = " ".join(str(pg.get(k, "")) for k in ("PN", "ID", "CL")).lower()
            if "bowtie2" in prog:
                return False
            if "bowtie" in prog:
                return True
    return False


def read_alignments(sam_path, xm_as_nh=False):
    """Iterates over primary mapped alignments. Excludes secondary (flag 256)."""
    mode = "rb" if str(sam_path).endswith(".bam") else "r"
    with pysam.AlignmentFile(sam_path, mode) as fh:
        for rec in fh:
            if _keep(rec):
                yield _to_aln(rec, xm_as_nh)


def read_region_handle(fh, chrom, start=None, end=None, xm_as_nh=False):
    """Iterate primary mapped alignments in a region of an ALREADY-OPEN AlignmentFile.

    Use this when fetching many regions from the same BAM: opening the file (which
    re-reads the header and index) once and reusing the handle avoids the per-region
    open/close overhead that dominates when there are tens of thousands of loci.
    """
    for rec in fh.fetch(chrom, start, end):
        if _keep(rec):
            yield _to_aln(rec, xm_as_nh)


def read_region(bam_path, chrom, start=None, end=None, xm_as_nh=False):
    """Iterates primary mapped alignments in a region of an INDEXED bam (fetch)."""
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        yield from read_region_handle(fh, chrom, start, end, xm_as_nh)


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


def reference_lengths(bam_path):
    """Reference name -> length, from the BAM header."""
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        return dict(zip(fh.references, fh.lengths))


def total_mapped(bam_path):
    """Total number of mapped alignments (fast, from the BAM index)."""
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        return fh.mapped


def weighted_library_size(sam_path, xm_as_nh=False):
    """Library size on the 1/NH scale: sum of the per-alignment weights.

    This is the number of distinct reads, not of alignment lines, so RPM lands on
    the same scale as the per-locus counts. Costs a full pass over the file (the
    index only knows alignment lines), unlike `total_mapped`.
    """
    return sum(1.0 / a.nh for a in read_alignments(sam_path, xm_as_nh))


def reference_alignment_counts(bam_path):
    """Alignments per reference (from the index, no scan).

    Sets the memory footprint of `classify`: each worker loads one whole
    reference at a time, so N workers hold the N references at once.
    """
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        return [s.mapped for s in fh.get_index_statistics() if s.mapped > 0]


def mapped_genome_size(bam_path):
    """Total length of the references that carry reads (from the index)."""
    with pysam.AlignmentFile(bam_path, "rb") as fh:
        return sum(fh.get_reference_length(s.contig)
                   for s in fh.get_index_statistics() if s.mapped > 0)
