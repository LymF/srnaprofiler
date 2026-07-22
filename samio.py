"""Reading of alignments (SAM/BAM) for the classification engine."""
from collections import namedtuple
import pysam

Aln = namedtuple("Aln", ["chrom", "pos", "strand", "seq", "nh"])


def read_alignments(sam_path):
    """Iterates over primary mapped alignments. Excludes secondary (flag 256)."""
    mode = "rb" if str(sam_path).endswith(".bam") else "r"
    with pysam.AlignmentFile(sam_path, mode) as fh:
        for rec in fh:
            if rec.is_unmapped or rec.is_secondary or rec.is_supplementary:
                continue
            seq = rec.query_sequence or ""
            if not seq:
                continue
            nh = rec.get_tag("NH") if rec.has_tag("NH") else 1
            strand = "-" if rec.is_reverse else "+"
            yield Aln(rec.reference_name, rec.reference_start, strand, seq, nh)
