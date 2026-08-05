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


def extend_and_merge(islands, pad, chrom_len=None):
    """Grow each island by `pad` on both sides, then merge what overlaps.

    Effective tolerated gap is 2*pad.
    """
    if not islands:
        return []
    grown = []
    for s, e in sorted(islands):
        s = max(0, s - pad)
        e = e + pad if chrom_len is None else min(chrom_len, e + pad)
        grown.append((s, e))
    merged = [grown[0]]
    for s, e in grown[1:]:
        ls, le = merged[-1]
        if s <= le:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged


def _compatible(a, b, max_strand_diff):
    """Whether two neighbouring profiles may be merged."""
    if a["size"] is not None and b["size"] is not None and a["size"] != b["size"]:
        return False
    ta, tb = a["frac_top"], b["frac_top"]
    if ta is not None and tb is not None and abs(ta - tb) > max_strand_diff:
        return False
    return True


def merge_compatible(islands, profiles, pad, chrom_len=None, max_strand_diff=0.5):
    """Extend and merge, but only across neighbours with a compatible profile.

    Unknown values (None) never block a merge.
    """
    if not islands:
        return []
    order = sorted(range(len(islands)), key=lambda i: islands[i])
    out = []
    cur_s, cur_e = None, None
    cur_prof = None
    for i in order:
        s, e = islands[i]
        s = max(0, s - pad)
        e = e + pad if chrom_len is None else min(chrom_len, e + pad)
        prof = profiles[i]
        if cur_s is None:
            cur_s, cur_e, cur_prof = s, e, prof
            continue
        if s <= cur_e and _compatible(cur_prof, prof, max_strand_diff):
            cur_e = max(cur_e, e)
            if cur_prof["size"] is None:
                cur_prof = prof
        else:
            out.append((cur_s, cur_e))
            cur_s, cur_e, cur_prof = s, e, prof
    out.append((cur_s, cur_e))
    return out


def cluster_regions(sam_path, mincov, pad):
    """Builds per-chromosome coverage and calls islands.

    Same rule as the per-chromosome path in `features`: raw islands, then a
    profile-conditioned bilateral merge.
    """
    from features import island_profiles

    depth = defaultdict(lambda: defaultdict(int))
    by_chrom = defaultdict(list)
    for a in read_alignments(sam_path):
        by_chrom[a.chrom].append(a)
        for p in range(a.pos, a.pos + len(a.seq)):
            depth[a.chrom][p] += 1
    regions = []
    for chrom in sorted(depth):
        islands = call_islands(depth[chrom], mincov, 0)
        profiles = island_profiles(islands, by_chrom[chrom])
        for start, end in merge_compatible(islands, profiles, pad):
            regions.append((chrom, start, end))
    return regions
