"""Role-separated markdown report built from per-library metrics.json files."""
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from benchmark.manifest import Library

DASH = "—"
STATUS_OK, STATUS_STALE, STATUS_INVALID = "ok", "stale", "INVALID"


@dataclass
class LibraryResult:
    lib: Library
    metrics: dict
    meta: dict


def collect_results(results_dir: Path, libs: list[Library]) -> list[LibraryResult]:
    """Load every <results_dir>/<ACC>/metrics.json whose accession is in libs."""
    by_acc = {l.accession: l for l in libs}
    out = []
    for mpath in sorted(Path(results_dir).glob("*/metrics.json")):
        metrics = json.loads(mpath.read_text())
        lib = by_acc.get(metrics.get("accession", mpath.parent.name))
        meta_path = mpath.parent / "benchmark_run.json"
        if lib is None or not meta_path.exists():
            continue
        out.append(LibraryResult(lib, metrics, json.loads(meta_path.read_text())))
    return out


def _status(r: LibraryResult, rules_dir: Path) -> str:
    rules = Path(rules_dir) / Path(r.meta["rules_path"]).name
    current = hashlib.sha256(rules.read_bytes()).hexdigest() if rules.exists() else None
    if current == r.meta.get("rules_sha256"):
        return STATUS_OK
    return STATUS_INVALID if r.lib.role == "test" else STATUS_STALE


def _fmt(v, digits: int = 3) -> str:
    if v is None:
        return DASH
    if isinstance(v, float):
        return f"{v:.{digits}f}"
    return str(v)


def _table(header: list[str], rows: list[list[str]]) -> list[str]:
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return lines + [""]


def _num_sorted(d: dict) -> list:
    return sorted(d, key=lambda k: int(k))


def _classes(items: list[dict]) -> list[str]:
    return sorted({c for per_class in items for c in per_class})


def _class_cells(per_class: dict, classes: list[str]) -> list[str]:
    cells = []
    for c in classes:
        m = per_class.get(c, {})
        cells += [_fmt(m.get("precision")), _fmt(m.get("recall")), _fmt(m.get("fdr"))]
    return cells


def _class_header(classes: list[str]) -> list[str]:
    return [f"{c} {k}" for c in classes for k in ("precision", "recall", "fdr")]


def _species_section(species: str, items: list[tuple[LibraryResult, str]]) -> list[str]:
    lines = [f"### {species}", ""]
    per_class = [r.metrics["classification"]["per_class"] for r, _ in items]
    classes = _classes(per_class)
    lines += ["Classification", ""]
    lines += _table(
        ["library", "status"] + _class_header(classes) + ["ambiguous_coverage", "ece", "unknown"],
        [[r.lib.accession, st] + _class_cells(r.metrics["classification"]["per_class"], classes)
         + [_fmt(r.metrics["classification"]["ambiguous_coverage"]),
            _fmt(r.metrics["classification"]["ece"]),
            _fmt(r.metrics["classification"]["unknown"]["count"])]
         for r, st in items])

    tiny = _num_sorted(items[0][0].metrics["delimitation"]["tiny_counts"])
    low = _num_sorted(items[0][0].metrics["delimitation"]["low_read_counts"])
    lines += ["Delimitation", ""]
    rows = []
    for r, st in items:
        d = r.metrics["delimitation"]
        s = d["summary"]
        rows.append([r.lib.accession, st, _fmt(s["n_clusters"]),
                     _fmt(s["mean_recovered_fraction"]), _fmt(s["mean_fragments"]),
                     _fmt(d["merged_distinct"])]
                    + [_fmt(d["tiny_counts"].get(k)) for k in tiny]
                    + [_fmt(d["low_read_counts"].get(k)) for k in low])
    lines += _table(["library", "status", "n_clusters", "mean_recovered_fraction",
                     "mean_fragments", "merged_distinct"]
                    + [f"tiny<={k}" for k in tiny] + [f"reads<{k}" for k in low], rows)

    sources = sorted({s for r, _ in items for s in r.metrics["by_source"]})
    for src in sources:
        with_src = [(r, st) for r, st in items if src in r.metrics["by_source"]]
        classes = _classes([r.metrics["by_source"][src] for r, _ in with_src])
        lines += [f"Source: {src}", ""]
        lines += _table(["library", "status"] + _class_header(classes),
                        [[r.lib.accession, st] + _class_cells(r.metrics["by_source"][src], classes)
                         for r, st in with_src])
    return lines


def _role_section(title: str, items: list[tuple[LibraryResult, str]]) -> list[str]:
    lines = [f"## {title}", ""]
    if not items:
        return lines + ["No results.", ""]
    for species in sorted({r.lib.organism for r, _ in items}):
        lines += _species_section(species, [(r, st) for r, st in items if r.lib.organism == species])
    return lines


def _aggregate(valid: list[LibraryResult]) -> list[str]:
    lines = ["### Aggregate", ""]
    if not valid:
        return lines + ["no valid test libraries", ""]
    lines += [f"Per-class mean over {len(valid)} valid test libraries (libraries without a value skipped).", ""]
    classes = _classes([r.metrics["classification"]["per_class"] for r in valid])
    rows = []
    for c in classes:
        cells = []
        for k in ("precision", "recall", "fdr"):
            vals = [r.metrics["classification"]["per_class"].get(c, {}).get(k) for r in valid]
            vals = [v for v in vals if v is not None]
            cells.append(_fmt(sum(vals) / len(vals)) if vals else DASH)
        rows.append([c] + cells)
    return lines + _table(["class", "precision", "recall", "fdr"], rows)


def build_report(results: list[LibraryResult], rules_dir: Path) -> str:
    """Markdown report; tuning and test results are never merged."""
    statused = [(r, _status(r, rules_dir)) for r in sorted(results, key=lambda r: r.lib.accession)]
    tuning = [(r, st) for r, st in statused if r.lib.role == "tuning"]
    test = [(r, st) for r, st in statused if r.lib.role == "test"]
    lines = ["# sRNAProfiler benchmark report", ""]
    lines += _role_section("Ajuste", tuning)
    lines += _role_section("Teste", test)
    lines += _aggregate([r for r, st in test if st == STATUS_OK])
    return "\n".join(lines)
