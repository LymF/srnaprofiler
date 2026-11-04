"""Command line entry point: python -m benchmark <command>."""
import argparse
import sys
from pathlib import Path

from benchmark.manifest import load_manifest, validate_coverage

DEFAULT_MANIFEST = Path(__file__).resolve().parent / "manifest.tsv"


def check_manifest(args: argparse.Namespace) -> int:
    try:
        libs = load_manifest(args.manifest)
        validate_coverage(libs)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"OK: {len(libs)} libraries, coverage satisfied")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="benchmark")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("check-manifest", help="load and validate the manifest")
    p.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    p.set_defaults(func=check_manifest)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
