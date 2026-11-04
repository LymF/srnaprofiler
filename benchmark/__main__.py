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


def report(args: argparse.Namespace) -> int:
    from benchmark.report import build_report, collect_results
    results = collect_results(args.results_dir, load_manifest(args.manifest))
    text = build_report(results, REPO)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(text)
    print(args.out)
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
    p = sub.add_parser("report", help="build the role-separated markdown report")
    p.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.add_argument("--out", type=Path, default=DEFAULT_RESULTS / "REPORT.md")
    p.set_defaults(func=report)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
