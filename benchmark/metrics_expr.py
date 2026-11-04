"""Expressed truth: which truth regions have enough reads in the library BAM."""
from pathlib import Path

import pysam

from benchmark.truth import Region


def expressed_regions(bam: Path, regions: list[Region], min_reads: int) -> set[int]:
    """Return id(region) for each region with >= min_reads alignments overlapping it.

    Counting choice: every alignment record stored in the BAM that overlaps the region is
    counted (pysam `count`, no filter), so primary and secondary alignments both count, and
    a multi-mapping read adds one per reported alignment. Unmapped reads are not in the
    index and are not counted. Contigs absent from the BAM count 0.
    """
    expressed: set[int] = set()
    with pysam.AlignmentFile(str(bam), "rb") as af:
        known = set(af.references)
        for r in regions:
            if r.chrom in known and af.count(r.chrom, r.start, r.end, read_callback="nofilter") >= min_reads:
                expressed.add(id(r))
    return expressed
