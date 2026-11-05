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


def profile_from_mass(mass, plus, window_mass, size, min_reads=0,
                      min_window_frac=0.5):
    """Merge-rule profile of a read mass, for one island or a running cluster.

    One definition for both, so a cluster of islands reads like one island of
    the same reads. Below `min_reads` (or empty) the profile is unknown, and
    unknown never blocks a merge. `window` is the size window holding more than
    `min_window_frac` of the mass; the mass keys travel along so profiles can be
    added up (`_add_profiles`).
    """
    prof = {"size": None, "window": None, "frac_top": None,
            "mass": mass, "plus": plus, "window_mass": window_mass}
    if not mass or mass < min_reads:
        return prof
    prof["size"] = size
    prof["frac_top"] = round(plus / mass, 4)
    if window_mass:
        best = max(range(len(window_mass)), key=lambda i: window_mass[i])
        if window_mass[best] / mass > min_window_frac:
            prof["window"] = best
    return prof


def _add_profiles(a, b, min_reads=0, min_window_frac=0.5):
    """Profile of the read mass of `a` and `b` together.

    Islands below `min_reads` have an unknown profile but their reads are still
    reads of the locus, so they count here. The size is the first one known
    (the cluster's first island with a call): it is only used when no window
    applies, and the modal length of a merged mass is no steadier than a
    single island's.
    """
    wa, wb = a.get("window_mass"), b.get("window_mass")
    wm = [x + y for x, y in zip(wa, wb)] if wa and wb else (wa or wb)
    size = a["size"] if a["size"] is not None else b["size"]
    return profile_from_mass(a.get("mass", 0.0) + b.get("mass", 0.0),
                             a.get("plus", 0.0) + b.get("plus", 0.0), wm, size,
                             min_reads, min_window_frac)


def _opposite_strands(ta, tb, cutoff):
    """Both fractions are stranded (ShortStack's call) and on opposite strands.

    ShortStack calls '+' at FracTop >= cutoff, '-' at <= 1 - cutoff, '.' between.
    An unstranded side ('.') is no evidence of a different strand.
    """
    # eps: 1 - 0.8 is 0.19999999999999996, which would leave a fraction of
    # exactly 0.2 unstranded.
    lo = 1 - cutoff + 1e-9
    return (ta >= cutoff and tb <= lo) or (tb >= cutoff and ta <= lo)


def _compatible(a, b, strand_cutoff, size_windows=None, strand_min_reads=None):
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
    if (ta is not None and tb is not None
            and _opposite_strands(ta, tb, strand_cutoff)
            and strand_min_reads is not None
            # Few reads make a strand fraction noise, not a strand: it blocks
            # only when both sides have enough mass to mean it.
            and a.get("mass", 0.0) >= strand_min_reads
            and b.get("mass", 0.0) >= strand_min_reads):
        return False
    return True


def merge_compatible(islands, profiles, pad, chrom_len=None, strand_cutoff=0.8,
                     size_windows=None, strand_min_reads=None,
                     profile_min_reads=0, min_window_frac=0.5):
    """Extend and merge, but only across neighbours with a compatible profile.

    Unknown values (None) never block a merge. A new island is compared with
    the profile of everything merged into the running cluster so far, not with
    its first island: a locus drifts, and the first island says little about
    its far end. The strand test blocks only two sides stranded at `strand_cutoff`
    on opposite strands; `strand_min_reads` None means it never blocks.
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
        if s <= cur_e and _compatible(cur_prof, prof, strand_cutoff,
                                      size_windows, strand_min_reads):
            cur_e, cur_raw_e = max(cur_e, e), max(cur_raw_e, raw_e)
            cur_prof = _add_profiles(cur_prof, prof, profile_min_reads,
                                     min_window_frac)
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


def carve_out(regions, intervals, min_piece=20):
    """Split regions so each known interval becomes a locus of its own.

    A miRNA hairpin swallowed by a cluster is a small fraction of it, so the
    cluster's class is decided by everything around the hairpin. Carving keeps
    the flanks as their own loci instead of discarding them. Pieces shorter
    than `min_piece` are dropped; the carved interval itself always survives.
    """
    if not intervals:
        return list(regions)
    cuts = sorted(intervals)
    out = []
    for s, e in sorted(regions):
        pos = s
        for cs, ce in cuts:
            cs, ce = max(cs, s), min(ce, e)
            if cs >= ce or cs < pos:
                continue
            if cs > pos and cs - pos >= min_piece:
                out.append((pos, cs))
            out.append((cs, ce))
            pos = ce
        if pos < e and (e - pos >= min_piece or pos == s):
            out.append((pos, e))
    return out


def cluster_regions(sam_path, mincov, pad, strand_min_reads=None,
                    profile_min_reads=0, strand_cutoff=0.8):
    """Builds per-chromosome coverage and calls islands.

    Same rule as the per-chromosome path in `features`: raw islands, then a
    profile-conditioned bilateral merge (`strand_min_reads` and
    `profile_min_reads` as there; the strand test does not block by default).
    """
    from features import island_profiles

    depth = defaultdict(lambda: defaultdict(float))
    by_chrom = defaultdict(list)
    for a in read_alignments(sam_path):
        by_chrom[a.chrom].append(a)
        w = 1.0 / a.nh
        for p in range(a.pos, a.pos + len(a.seq)):
            depth[a.chrom][p] += w
    regions = []
    for chrom in sorted(depth):
        islands = call_islands(depth[chrom], mincov, 0)
        profiles = island_profiles(islands, by_chrom[chrom],
                                   min_reads=profile_min_reads)
        for start, end in merge_compatible(
                islands, profiles, pad, strand_min_reads=strand_min_reads,
                strand_cutoff=strand_cutoff,
                profile_min_reads=profile_min_reads):
            regions.append((chrom, start, end))
    return regions
