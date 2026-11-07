#!/usr/bin/env python3
"""Build negative-control regions (rRNA and protein-coding exons) for one genome.

Usage: build_controls.py ANNOT.gff3[.gz] TRUTH [TRUTH ...] [--aliases TSV] -o controls.bed

Controls are the union of
  * features of type `rRNA` (and `rRNA_gene`), written with name "rRNA";
  * `exon` features whose Parent is an `mRNA` (Ensembl and NCBI both label protein-coding
    transcripts mRNA), written with name "exon";
minus every control interval that overlaps (>0 bp) any region of any TRUTH file. TRUTH files
are BED (0-based) or GFF3 (.gff/.gff3, optionally .gz; 1-based), any class.
--aliases is a TSV: column 1 = contig name in the annotation, column 2 = genome contig name;
annotation contigs without an alias keep their name. Output: 4-column BED, sorted.
"""
import argparse
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from benchmark.metrics_class import IntervalIndex  # noqa: E402

RRNA_TYPES = {"rRNA", "rRNA_gene"}


def _open(path):
    p = str(path)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def _attrs(field: str) -> dict:
    out = {}
    for kv in field.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k] = v
    return out


def _is_gff(path) -> bool:
    return str(path).removesuffix(".gz").lower().endswith((".gff", ".gff3"))


def read_truth(path) -> list[tuple[str, int, int]]:
    out = []
    gff = _is_gff(path)
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if gff:
                out.append((f[0], int(f[3]) - 1, int(f[4])))
            else:
                out.append((f[0], int(f[1]), int(f[2])))
    return out


def read_aliases(path) -> dict[str, str]:
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] and not line.startswith("#"):
                out[f[0]] = f[1]
    return out


def read_controls(gff, aliases: dict[str, str]) -> list[tuple[str, int, int, str]]:
    """Deduplicated (chrom, start0, end, name) control candidates from a GFF3."""
    mrna, exons, out = set(), [], set()
    with _open(gff) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom = aliases.get(f[0], f[0])
            typ = f[2]
            a = _attrs(f[8])
            if typ == "mRNA" and "ID" in a:
                mrna.add(a["ID"])
            elif typ == "exon":
                exons.append((chrom, int(f[3]) - 1, int(f[4]), a.get("Parent", "").split(",")))
            elif typ in RRNA_TYPES:
                out.add((chrom, int(f[3]) - 1, int(f[4]), "rRNA"))
    for chrom, s, e, parents in exons:
        if any(p in mrna for p in parents):
            out.add((chrom, s, e, "exon"))
    return sorted(out)


def build(gff, truths, aliases=None) -> list[tuple[str, int, int, str]]:
    cands = read_controls(gff, aliases or {})
    truth = [iv for t in truths for iv in read_truth(t)]
    if not truth:
        return cands
    index = IntervalIndex(truth)
    return [c for c in cands if next(index.overlapping(c[0], c[1], c[2]), None) is None]


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("gff3", help="annotation GFF3 (optionally .gz)")
    ap.add_argument("truth", nargs="+", help="truth regions to exclude (BED or GFF3)")
    ap.add_argument("--aliases", help="TSV: annotation contig -> genome contig")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args(argv)
    aliases = read_aliases(args.aliases) if args.aliases else {}
    rows = build(args.gff3, args.truth, aliases)
    with open(args.output, "w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"{len(rows)} control intervals -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
