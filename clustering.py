"""Cluster definition via coverage islands (ShortStack-style)."""
from collections import defaultdict
from samio import read_alignments


def call_islands(depth_by_pos, mincov, pad):
    """Positions with depth>=mincov, merging islands whose gap (empty span) <= pad.

    Returns [(start, end)] with end exclusive.
    """
    positions = sorted(p for p, d in depth_by_pos.items() if d >= mincov)
    if not positions:
        return []
    islands = []
    start = prev = positions[0]
    for p in positions[1:]:
        if p - prev - 1 > pad:
            islands.append((start, prev + 1))
            start = p
        prev = p
    islands.append((start, prev + 1))
    return islands


def cluster_regions(sam_path, mincov, pad):
    """Builds per-chromosome coverage and calls islands."""
    depth = defaultdict(lambda: defaultdict(int))
    for a in read_alignments(sam_path):
        for p in range(a.pos, a.pos + len(a.seq)):
            depth[a.chrom][p] += 1
    regions = []
    for chrom in sorted(depth):
        for start, end in call_islands(depth[chrom], mincov, pad):
            regions.append((chrom, start, end))
    return regions
