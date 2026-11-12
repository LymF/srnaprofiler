"""Data layer for the HTML report: read clusters.tsv and apply the depth floor."""
import csv
import functools
import os
import statistics
from dataclasses import dataclass, field

import yaml

_E_KEYS = ("1u", "5a", "fractop_bid", "pingpong", "phase")
# Classes only the plant config has; their gate_* columns or calls mark a plant run.
PLANT_CLASSES = ("phasiRNA", "phasiRNA24", "hcsiRNA")
_RULES = {"animal": "rules.yaml", "plant": "rules.plant.yaml"}
# defining_evidence names in the rules files -> report check keys
_EV_CHECK = {"1u": "1u", "5a": "5a", "fractop_bid": "strand", "pingpong": "pingpong",
             "phase": "phasing", "size": "size"}


def _num(x):
    """float or None (empty, NA, N and unparsable values)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _int(x):
    v = _num(x)
    return None if v is None else int(v)


def _str(x):
    return "" if x is None or x in ("NA", "N") else x


@dataclass
class Cluster:
    name: str
    chrom: str
    start: int
    end: int
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
    pct_5a: float = None


def _row_to_cluster(row, has_evidence):
    n = _num(row.get("dicer_call"))
    return Cluster(
        name=row.get("name", ""), chrom=row.get("chrom", ""),
        start=_int(row.get("start")), end=_int(row.get("end")),
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
        pct_5a=_num(row.get("pct_5a")),
    )


def load_clusters(path):
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        has_ev = all(f"e_{k}" in (rd.fieldnames or []) for k in _E_KEYS)
        return [_row_to_cluster(r, has_ev) for r in rd]


def detect_kingdom(clusters):
    """"plant" when any plant-only class shows up as a gate_* column or a call.

    The TSV does not record which rules file classify used; the plant classes are
    the only reliable trace of it. An old TSV with neither reads as "animal".
    """
    for c in clusters:
        if c.cls in PLANT_CLASSES or any(k in PLANT_CLASSES for k in (c.gates or {})):
            return "plant"
    return "animal"


def detect_mode(clusters):
    """"ml" when classify ran with --classifier ml (evidence column "ml:<algorithm>")."""
    return "ml" if any((c.evidence or "").startswith("ml:") for c in clusters) else "rules"


@functools.lru_cache(maxsize=None)
def class_specs(kingdom):
    """Per-class display facts from the shipped rules file of this kingdom.

    {class: {"defining": set of check keys, "size": (lo, hi) or None,
    "min_length": bp or None, "phase_max_p": float or None, "by_phase": bool}}.
    A run with a custom --config may differ; the report cannot know it.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), _RULES[kingdom])
    with open(path) as fh:
        classes = (yaml.safe_load(fh) or {}).get("classes") or {}
    out = {}
    for name, spec in classes.items():
        spec = spec or {}
        if not spec.get("enabled", True):
            continue
        names = set()
        for slot in spec.get("defining_evidence", []):
            names.update(slot if isinstance(slot, list) else [slot])
        size = spec.get("size") or {}
        out[name] = {"defining": {_EV_CHECK[n] for n in names if n in _EV_CHECK},
                     "size": (size["lo"], size["hi"]) if "lo" in size and "hi" in size else None,
                     "min_length": spec.get("cluster_min_length"),
                     "phase_max_p": spec.get("phase_max_pvalue"),
                     "by_phase": "phase" in names}
    return out


@dataclass
class Selection:
    shown: list
    detail: list
    n_total: int
    n_hidden: int
    floor: float
    capped: bool
    kingdom: str = "animal"
    mode: str = "rules"


def select(clusters, floor, max_detail=3000):
    shown = [c for c in clusters if c.reads is not None and c.reads >= floor]
    # ambiguous rows carry the nearest class's confidence: list them after resolved ones
    shown.sort(key=lambda c: (c.cls == "ambiguous",
                              -(c.confidence if c.confidence is not None else -1.0), -c.reads))
    return Selection(shown=shown, detail=shown[:max_detail], n_total=len(clusters),
                     n_hidden=len(clusters) - len(shown), floor=floor,
                     capped=len(shown) > max_detail, kingdom=detect_kingdom(clusters),
                     mode=detect_mode(clusters))


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


def _percentiles(values):
    """0-100 int percentile rank per value; ties share the mean rank. One value: 100."""
    n = len(values)
    if n == 1:
        return [100]
    order = sorted(values)
    first, last = {}, {}
    for i, v in enumerate(order):
        first.setdefault(v, i)
        last[v] = i
    return [round(100 * (first[v] + last[v]) / 2 / (n - 1)) for v in values]


def library_context(clusters):
    """{name: {"reads_pct", "conf_pct"}}: where a cluster sits among the shown ones.

    reads_pct ranks reads over all clusters; conf_pct ranks confidence within the
    cluster's class, None for ambiguous rows, a missing confidence or a class of one.
    """
    ok = [c for c in clusters if c.reads is not None]
    out = {c.name: {"reads_pct": None, "conf_pct": None} for c in clusters}
    if ok:
        for c, p in zip(ok, _percentiles([c.reads for c in ok])):
            out[c.name]["reads_pct"] = p
    by_cls = {}
    for c in clusters:
        if c.cls != "ambiguous" and c.confidence is not None:
            by_cls.setdefault(c.cls, []).append(c)
    for grp in by_cls.values():
        if len(grp) > 1:
            for c, p in zip(grp, _percentiles([c.confidence for c in grp])):
                out[c.name]["conf_pct"] = p
    return out
