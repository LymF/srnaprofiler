"""Role-separated markdown report built from per-library metrics.json files."""
import json
from dataclasses import dataclass
from pathlib import Path

from benchmark.manifest import Library
from benchmark.run import pipeline_sha256

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


def _status(r: LibraryResult, repo_dir: Path) -> str:
    """ok only if the result was made by the pipeline code and rules in repo_dir now."""
    if r.meta.get("pipeline_sha256") == pipeline_sha256(repo_dir):
        return STATUS_OK
    return STATUS_INVALID if r.lib.role == "test" else STATUS_STALE


def _short_sha(r: LibraryResult) -> str:
    sha = r.meta.get("pipeline_sha256")
    return sha[:8] if sha else DASH


def _mirna_cell(r: LibraryResult) -> str:
    used = r.meta.get("mirna_loci_used")
    return DASH if used is None else ("yes" if used else "no")


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


CLASS_COLS = [("precision", "precision (annotated loci)"), ("tp", "tp"), ("fp", "fp"),
              ("n_called", "n_called"), ("n_unknown", "n_unknown"), ("recall", "recall"),
              ("recall_expressed", "recall_expressed"), ("n_truth", "n_truth"),
              ("n_truth_expressed", "n_truth_expressed"), ("fdr", "fdr")]


def _class_cells(per_class: dict, classes: list[str]) -> list[str]:
    return [_fmt(per_class.get(c, {}).get(key)) for c in classes for key, _ in CLASS_COLS]


def _class_header(classes: list[str]) -> list[str]:
    return [f"{c} {label}" for c in classes for _, label in CLASS_COLS]


def _species_section(species: str, items: list[tuple[LibraryResult, str]]) -> list[str]:
    lines = [f"### {species}", ""]
    per_class = [r.metrics["classification"]["per_class"] for r, _ in items]
    classes = _classes(per_class)
    lines += ["Classification", ""]
    lines += _table(
        ["library", "status", "pipeline", "mirna_loci_used"] + _class_header(classes)
        + ["ambiguous_coverage", "ece", "unknown"],
        [[r.lib.accession, st, _short_sha(r), _mirna_cell(r)] + _class_cells(r.metrics["classification"]["per_class"], classes)
         + [_fmt(r.metrics["classification"]["ambiguous_coverage"]),
            _fmt(r.metrics["classification"]["ece"]),
            _fmt(r.metrics["classification"]["unknown"]["count"])]
         for r, st in items])

    tiny = _num_sorted(items[0][0].metrics["delimitation"]["global"]["tiny_counts"])
    low = _num_sorted(items[0][0].metrics["delimitation"]["global"]["low_read_counts"])
    lines += ["Loci counts (all loci, independent of truth source)", ""]
    lines += _table(["library", "status"] + [f"tiny<={k}" for k in tiny] + [f"reads<{k}" for k in low]
                    + ["fragments", "fragment_reads"],
                    [[r.lib.accession, st]
                     + [_fmt(r.metrics["delimitation"]["global"]["tiny_counts"].get(k)) for k in tiny]
                     + [_fmt(r.metrics["delimitation"]["global"]["low_read_counts"].get(k)) for k in low]
                     + [_fmt(r.metrics["delimitation"]["global"].get("n_fragments")),
                        _fmt(r.metrics["delimitation"]["global"].get("fragment_read_frac"))]
                     for r, st in items])
    delim_sources = sorted({s for r, _ in items for s in r.metrics["delimitation"]["by_source"]})
    for src in delim_sources:
        lines += [f"Delimitation: {src}", ""]
        rows = []
        for r, st in items:
            d = r.metrics["delimitation"]["by_source"].get(src)
            if d is None:
                continue
            s = d["summary"]
            rows.append([r.lib.accession, st, _fmt(s["n_clusters"]),
                         _fmt(s["mean_recovered_fraction"]), _fmt(s["mean_fragments"]),
                         _fmt(d["merged_distinct"])])
        lines += _table(["library", "status", "n_clusters", "mean_recovered_fraction",
                         "mean_fragments", "merged_distinct"], rows)

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


AGG_KEYS = ("precision", "recall", "recall_expressed", "fdr")


def _mean(vals: list) -> float | None:
    vals = [v for v in vals if v is not None]
    return sum(vals) / len(vals) if vals else None


def _aggregate(valid: list[LibraryResult]) -> list[str]:
    """Mean per class over libraries within a species, then over species within a taxon group.

    Arthropods and plants are never mixed, and every species weighs the same however many
    libraries it has.
    """
    lines = ["### Aggregate", ""]
    if not valid:
        return lines + ["no valid test libraries", ""]
    lines += [f"{len(valid)} valid test librar{'y' if len(valid) == 1 else 'ies'}. Per-class mean over "
              "the libraries of each species, then mean over the species of each taxon group "
              "(values missing at either level are skipped).", ""]
    for group in sorted({r.lib.taxon_group for r in valid}):
        in_group = [r for r in valid if r.lib.taxon_group == group]
        species = sorted({r.lib.organism for r in in_group})
        lines += [f"#### {group}", ""]
        lines += _table(["species", "libraries"],
                        [[sp, str(sum(1 for r in in_group if r.lib.organism == sp))] for sp in species])
        classes = _classes([r.metrics["classification"]["per_class"] for r in in_group])
        rows = []
        for c in classes:
            cells = []
            for k in AGG_KEYS:
                per_species = [
                    _mean([r.metrics["classification"]["per_class"].get(c, {}).get(k)
                           for r in in_group if r.lib.organism == sp])
                    for sp in species]
                cells.append(_fmt(_mean(per_species)))
            rows.append([c] + cells)
        lines += _table(["class", "precision (annotated loci)", "recall", "recall_expressed", "fdr"],
                        rows)
    return lines


def build_report(results: list[LibraryResult], repo_dir: Path) -> str:
    """Markdown report; tuning and test results are never merged."""
    statused = [(r, _status(r, repo_dir)) for r in sorted(results, key=lambda r: r.lib.accession)]
    tuning = [(r, st) for r, st in statused if r.lib.role == "tuning"]
    test = [(r, st) for r, st in statused if r.lib.role == "test"]
    lines = ["# sRNAProfiler benchmark report", ""]
    lines += _role_section("Ajuste", tuning)
    lines += _role_section("Teste", test)
    lines += _aggregate([r for r, st in test if st == STATUS_OK])
    return "\n".join(lines)
