"""Data layer for the HTML report: read clusters.tsv and apply the depth floor."""
import csv
import statistics
from dataclasses import dataclass, field

_E_KEYS = ("1u", "5a", "fractop_bid", "pingpong", "phase")


def _num(x):
    """float or None (empty, NA, N and unparsable values)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _str(x):
    return "" if x is None or x in ("NA", "N") else x


@dataclass
class Cluster:
    name: str
    chrom: str
    start: float
    end: float
    reads: float
    reads_rpm: float
    cls: str
    confidence: float
    dicer_call: int
    size_purity: float
    frac_top: float
    pct_1u: float
    pct_10a: float
    pingpong_z: float
    phase_pvalue: float
    multimap_frac: float
    unique_frac: float
    review: bool
    major_rna: str
    evidence: str
    e: dict = field(default_factory=dict)
    gates: dict = field(default_factory=dict)
    has_evidence: bool = False


def _row_to_cluster(row, has_evidence):
    n = _num(row.get("dicer_call"))
    return Cluster(
        name=row.get("name", ""), chrom=row.get("chrom", ""),
        start=_num(row.get("start")), end=_num(row.get("end")),
        reads=_num(row.get("reads")), reads_rpm=_num(row.get("reads_rpm")),
        cls=_str(row.get("cls")), confidence=_num(row.get("confidence")),
        dicer_call=None if n is None else int(n),
        size_purity=_num(row.get("size_purity")), frac_top=_num(row.get("frac_top")),
        pct_1u=_num(row.get("pct_1u")), pct_10a=_num(row.get("pct_10a")),
        pingpong_z=_num(row.get("pingpong_z")), phase_pvalue=_num(row.get("phase_pvalue")),
        multimap_frac=_num(row.get("multimap_frac")), unique_frac=_num(row.get("unique_frac")),
        review=row.get("review") == "True",
        major_rna=_str(row.get("major_rna")), evidence=_str(row.get("evidence")),
        e={k: _num(row.get(f"e_{k}")) for k in _E_KEYS},
        gates={k[5:]: v == "1" for k, v in row.items()
               if k.startswith("gate_") and _num(v) is not None},
        has_evidence=has_evidence,
    )


def load_clusters(path):
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        has_ev = all(f"e_{k}" in (rd.fieldnames or []) for k in _E_KEYS)
        return [_row_to_cluster(r, has_ev) for r in rd]


@dataclass
class Selection:
    shown: list
    detail: list
    n_total: int
    n_hidden: int
    floor: float
    capped: bool


def select(clusters, floor, max_detail=3000):
    shown = [c for c in clusters if c.reads is not None and c.reads >= floor]
    shown.sort(key=lambda c: (-(c.confidence if c.confidence is not None else -1.0), -c.reads))
    return Selection(shown=shown, detail=shown[:max_detail], n_total=len(clusters),
                     n_hidden=len(clusters) - len(shown), floor=floor,
                     capped=len(shown) > max_detail)


def library_summary(clusters, fragments_path):
    libs = [c.reads / c.reads_rpm * 1e6 for c in clusters
            if c.reads is not None and c.reads_rpm and c.reads_rpm > 0]
    wsum = sum(c.reads for c in clusters if c.reads is not None and c.unique_frac is not None)
    uniq = (sum(c.reads * c.unique_frac for c in clusters
                if c.reads is not None and c.unique_frac is not None) / wsum) if wsum else None
    n_frag = frag_reads = None
    if fragments_path:
        with open(fragments_path, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        n_frag = len(rows)
        frag_reads = sum(_num(r.get("reads")) or 0.0 for r in rows)
    return {"library_reads": statistics.median(libs) if libs else None,
            "cluster_reads": sum(c.reads for c in clusters if c.reads is not None),
            "unique_share": uniq, "fragments": n_frag, "fragment_reads": frag_reads}
