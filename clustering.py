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


def _same_size_group(sa, sb, size_windows):
    """Whether two dominant sizes belong to the same biogenesis window.

    Without windows, exact equality.
    """
    if not size_windows:
        return sa == sb
    ga = frozenset(i for i, (lo, hi) in enumerate(size_windows) if lo <= sa <= hi)
    gb = frozenset(i for i, (lo, hi) in enumerate(size_windows) if lo <= sb <= hi)
    # A size outside every window belongs to no known pathway, so it is not
    # evidence of a conflict -- same treatment as an unknown profile.
    if not ga or not gb:
        return True
    return bool(ga & gb)


def _compatible(a, b, max_strand_diff, size_windows=None):
    """Whether two neighbouring profiles may be merged."""
    wa, wb = a.get("window"), b.get("window")
    if wa is not None and wb is not None:
        # The window comes from the read mass, which holds up on the few reads a
        # raw island has; the modal length does not. Prefer it when both have one.
        if wa != wb:
            return False
    elif (wa is None and wb is None and a["size"] is not None
            and b["size"] is not None
            and not _same_size_group(a["size"], b["size"], size_windows)):
        return False
    ta, tb = a["frac_top"], b["frac_top"]
    if ta is not None and tb is not None and abs(ta - tb) > max_strand_diff:
        return False
    return True


def merge_compatible(islands, profiles, pad, chrom_len=None, max_strand_diff=0.5,
                     size_windows=None):
    """Extend and merge, but only across neighbours with a compatible profile.

    Unknown values (None) never block a merge.
    """
    if not islands:
        return []
    order = sorted(range(len(islands)), key=lambda i: islands[i])
    out = []
    cur_s = cur_e = cur_raw_e = None
    cur_prof = None
    for i in order:
        raw_s, raw_e = islands[i]
        s = max(0, raw_s - pad)
        e = raw_e + pad if chrom_len is None else min(chrom_len, raw_e + pad)
        prof = profiles[i]
        if cur_s is None:
            cur_s, cur_e, cur_raw_e, cur_prof = s, e, raw_e, prof
            continue
        if s <= cur_e and _compatible(cur_prof, prof, max_strand_diff, size_windows):
            cur_e, cur_raw_e = max(cur_e, e), max(cur_raw_e, raw_e)
            if cur_prof["size"] is None:
                cur_prof = prof
        else:
            # Not merged: they split the gap, because a read is assigned to
            # every locus it touches and overlapping loci would count it twice.
            if s < cur_e:
                bound = (cur_raw_e + raw_s) // 2
                cur_e, s = min(cur_e, bound), max(s, bound)
            out.append((cur_s, cur_e))
            cur_s, cur_e, cur_raw_e, cur_prof = s, e, raw_e, prof
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
