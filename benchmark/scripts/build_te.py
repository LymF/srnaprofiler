#!/usr/bin/env python3
"""Build the transposon BED used by the independent-biology metric from a UCSC rmsk table.

Usage: build_te.py rmsk.txt.gz -o rmsk_te.bed

Source (dm6): https://hgdownload.soe.ucsc.edu/goldenPath/dm6/database/rmsk.txt.gz
Keeps rows whose repClass is one of LTR, LINE, DNA, RC, Other, Unknown (exact match; uncertain
classes such as "LTR?" and non-TE classes such as Simple_repeat, Low_complexity, Satellite,
rRNA, tRNA are dropped). The `chr` prefix is stripped from contig names. Output: 3-column BED
(0-based start, as in rmsk), sorted. rmsk.txt columns used: genoName (6th), genoStart (7th),
genoEnd (8th), repClass (12th).
"""
import argparse
import gzip
import sys
from pathlib import Path

TE_CLASSES = {"LTR", "LINE", "DNA", "RC", "Other", "Unknown"}


def _open(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def build(rmsk) -> list[tuple[str, int, int]]:
    out = []
    with _open(rmsk) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12 or f[11] not in TE_CLASSES:
                continue
            chrom = f[5][3:] if f[5].startswith("chr") else f[5]
            out.append((chrom, int(f[6]), int(f[7])))
    return sorted(out)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("rmsk", help="UCSC rmsk.txt.gz")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args(argv)
    rows = build(args.rmsk)
    with open(args.output, "w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} TE intervals -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
