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


def delimitation_metrics(loci: list[dict], clusters: list[Region], cfg: Config) -> dict:
    """Score cluster delimitation; locus class is ignored (every locus counts)."""
    index = IntervalIndex([(l["chrom"], l["start"], l["end"]) for l in loci])
    regions_hit = [0] * len(loci)  # distinct truth regions each locus overlaps

    per_cluster = []
    for c in clusters:
        spans = []
        for idx, _ in index.overlapping(c.chrom, c.start, c.end):
            loc = loci[idx]
            spans.append((max(loc["start"], c.start), min(loc["end"], c.end)))
            regions_hit[idx] += 1
        per_cluster.append({
            "name": c.name, "chrom": c.chrom, "start": c.start, "end": c.end,
            "recovered_fraction": _union_bp(spans) / (c.end - c.start),
            "n_fragments": len(spans),
        })

    n = len(per_cluster)
    return {
        "per_cluster": per_cluster,
        "merged_distinct": sum(1 for h in regions_hit if h >= 2),
        "tiny_counts": {bp: sum(1 for l in loci if l["end"] - l["start"] <= bp)
                        for bp in cfg.tiny_bp},
        "low_read_counts": {r: sum(1 for l in loci if l["reads"] < r) for r in cfg.low_reads},
        "summary": {
            "n_clusters": n,
            "mean_recovered_fraction": sum(p["recovered_fraction"] for p in per_cluster) / n if n else 0.0,
            "mean_fragments": sum(p["n_fragments"] for p in per_cluster) / n if n else 0.0,
        },
    }
