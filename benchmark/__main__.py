"""Command line entry point: python -m benchmark <command>."""
import argparse
import os
import sys
from pathlib import Path

from benchmark.manifest import load_manifest, validate_coverage

REPO = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = Path(__file__).resolve().parent / "manifest.tsv"
DEFAULT_RESULTS = REPO / "benchmark" / "results"


def default_data_dir() -> Path:
    return Path(os.environ.get("BENCHMARK_DATA") or REPO / "data" / "benchmark")


def check_manifest(args: argparse.Namespace) -> int:
    try:
        libs = load_manifest(args.manifest)
        validate_coverage(libs)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"OK: {len(libs)} libraries, coverage satisfied")
    return 0


def run(args: argparse.Namespace) -> int:
    from benchmark.orchestrate import run_benchmark_library
    try:
        path = run_benchmark_library(args.library, args.manifest, args.data_dir,
                                     args.results_dir, args.threads)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return 1
    print(path)
    return 0


def metrics(args: argparse.Namespace) -> int:
    from benchmark.orchestrate import rebuild_metrics
    try:
        path = rebuild_metrics(args.library, args.manifest, args.results_dir)
    except (ValueError, FileNotFoundError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(path)
    return 0


def report(args: argparse.Namespace) -> int:
    from benchmark.report import build_report, collect_results
    results = collect_results(args.results_dir, load_manifest(args.manifest))
    text = build_report(results, REPO)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(args.out)
    return 0


def ml_cv(args: argparse.Namespace) -> int:
    from benchmark.ml.cv import cv_main
    try:
        n = cv_main(args.profile, [a for a in args.algorithms.split(",") if a],
                    [v for v in args.variants.split(",") if v], args.seeds,
                    args.results_dir, args.out, args.manifest, args.n_jobs)
    except (ValueError, FileNotFoundError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"{n} fold line(s) appended to {args.out / 'folds.jsonl'}")
    return 0


def ml_select(args: argparse.Namespace) -> int:
    from benchmark.ml.select import SelectionError, select
    algs = [a for a in args.algorithms.split(",") if a] if args.algorithms else None
    try:
        select(args.folds, args.output, algs, args.seeds)
    except SelectionError as exc:
        print(exc, file=sys.stderr)
        return 1
    print(args.output)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="benchmark")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check-manifest", help="load and validate the manifest")
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.set_defaults(func=check_manifest)
    p = sub.add_parser("run", help="run the pipeline and metrics for one library")
    p.add_argument("--library", required=True)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--data-dir", type=Path, default=default_data_dir())
    p.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    p.add_argument("--threads", type=int, default=8)
    p.set_defaults(func=run)
    p = sub.add_parser("metrics", help="rebuild metrics.json from an existing results/<ACC>/work")
    p.add_argument("--library", required=True)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    p.set_defaults(func=metrics)
    p = sub.add_parser("report", help="build the role-separated markdown report")
    p.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--out", type=Path, default=DEFAULT_RESULTS / "REPORT.md")
    p.set_defaults(func=report)
    p = sub.add_parser("ml-cv", help="cross-validate the ML classifiers (resumable)")
    p.add_argument("--profile", choices=["arthropod", "plant"], required=True)
    p.add_argument("--algorithms", default="all", help="comma list or 'all'")
    p.add_argument("--variants", default="main,stack,repro")
    p.add_argument("--seeds", type=int, default=5)
    p.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--n-jobs", type=int, default=8)
    p.set_defaults(func=ml_cv)
    p = sub.add_parser("ml-select", help="apply the selection rule to folds.jsonl")
    p.add_argument("--folds", type=Path, required=True)
    p.add_argument("-o", "--output", type=Path, required=True)
    p.add_argument("--seeds", type=int, default=None,
                   help="seeds the run must contain (default: the largest recorded in folds.jsonl)")
    p.add_argument("--algorithms", default=None, help="comma list to rank (default: all in file)")
    p.set_defaults(func=ml_select)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
