"""Delimitation metrics: how well predicted loci recover the extent of truth clusters."""
from benchmark.manifest import Config
from benchmark.metrics_class import IntervalIndex
from benchmark.truth import Region


def _union_bp(spans: list[tuple[int, int]]) -> int:
    """Total bp covered by the union of (start, end) spans."""
    total, cur_end = 0, None
    for s, e in sorted(spans):
        if cur_end is None or s > cur_end:
            total += e - s
            cur_end = e
        elif e > cur_end:
            total += e - cur_end
            cur_end = e
    return total


def _span(locus: dict) -> int:
    """Read span of a locus, falling back to its window when the run did not record one."""
    span = locus.get("read_span")
    return locus["end"] - locus["start"] if span is None else span


def locus_counts(loci: list[dict], cfg: Config, fragments: list[dict] | None = None) -> dict:
    """Tiny and low-read locus counts over all loci (not tied to any truth source).

    Tiny counts use the read span and skip miRNA loci (hairpin loci are legitimately short).
    fragments is the parsed fragments.tsv (None for runs that predate it): n_fragments and
    fragment_read_frac are None then. fragment_read_frac = fragment read mass / (fragment + locus read mass), miRNA loci included (which dilutes it in miRNA-rich libraries); reads overlapping two regions count in both, so it is a share of bucketed mass, not of library mass.
    """
    n_frag = frac = None
    if fragments is not None:
        n_frag = len(fragments)
        frag_reads = sum(f["reads"] for f in fragments)
        total = frag_reads + sum(l["reads"] for l in loci)
        frac = frag_reads / total if total else None
    non_mirna = [l for l in loci if l.get("cls") != "miRNA"]
    return {
        "tiny_counts": {bp: sum(1 for l in non_mirna if _span(l) <= bp) for bp in cfg.tiny_bp},
        "low_read_counts": {r: sum(1 for l in loci if l["reads"] < r) for r in cfg.low_reads},
        "n_fragments": n_frag,
        "fragment_read_frac": frac,
    }


def _joins_distinct(regions: list[Region], gap_bp: int) -> bool:
    """True if two regions are more than gap_bp apart (overlapping ones are never distinct)."""
    ordered = sorted(regions, key=lambda r: r.start)
    min_end = ordered[0].end
    for r in ordered[1:]:
        if r.start - min_end > gap_bp:
            return True
        min_end = min(min_end, r.end)
    return False


def delimitation_metrics(loci: list[dict], clusters: list[Region], cfg: Config) -> dict:
    """Score cluster delimitation against the clusters of ONE truth source.

    Locus class is ignored (every locus counts). Pass one source at a time: clusters of
    different sources overlap each other and would be double counted. A locus counts in
    merged_distinct when it overlaps two clusters more than cfg.merge_gap_bp apart.
    """
    index = IntervalIndex([(l["chrom"], l["start"], l["end"]) for l in loci])
    regions_hit: list[list[Region]] = [[] for _ in loci]  # clusters each locus overlaps

    per_cluster = []
    for c in clusters:
        spans = []
        for idx, _ in index.overlapping(c.chrom, c.start, c.end):
            loc = loci[idx]
            spans.append((max(loc["start"], c.start), min(loc["end"], c.end)))
            regions_hit[idx].append(c)
        per_cluster.append({
            "name": c.name, "chrom": c.chrom, "start": c.start, "end": c.end,
            "recovered_fraction": _union_bp(spans) / (c.end - c.start),
            "n_fragments": len(spans),
        })

    n = len(per_cluster)
    return {
        "per_cluster": per_cluster,
        "merged_distinct": sum(1 for h in regions_hit
                               if len(h) >= 2 and _joins_distinct(h, cfg.merge_gap_bp)),
        "summary": {
            "n_clusters": n,
            "mean_recovered_fraction": sum(p["recovered_fraction"] for p in per_cluster) / n if n else 0.0,
            "mean_fragments": sum(p["n_fragments"] for p in per_cluster) / n if n else 0.0,
        },
    }
