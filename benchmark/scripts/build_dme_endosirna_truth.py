#!/usr/bin/env python3
"""Build benchmark/truth/dmel/dme_endosiRNA_2008_dm6.bed from Czech et al. 2008 (Nature 453:798) Tables S1/S2.

Usage (from the repo root, CrossMap and pdftotext on PATH):
  python benchmark/scripts/build_dme_endosirna_truth.py SUPPL_PDF dm3ToDm6.over.chain.gz OUT.bed [WORKDIR]

SUPPL_PDF is https://static-content.springer.com/esm/art%3A10.1038%2Fnature07007/MediaObjects/41586_2008_BFnature07007_MOESM268_ESM.pdf
The chain is https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz.
Steps: pdftotext -layout; parse the two tables (Release 5 = dm3 coordinates, size == stop - start so they are BED-like);
drop transgene ("vector") rows; union overlapping loci of the two tables; CrossMap bed dm3 -> dm6; rename chrN -> N
(Ensembl BDGP6 names); drop loci overlapping miRNA (MirGeneDB, miRBase) or piRNA-cluster truth already in the repo.
Counts go to stderr. The downloaded supplement stays outside git.
"""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from benchmark.truth import load_aliases, load_regions  # noqa: E402

ROW = re.compile(r"^\s*(X|2L|2R|3L|3R|4)\s+(\d+)\s+(\d+)\s+(.+?)\s{2,}(\d+)\s+[\d.]+\s+\d+\s+\d+\s+\d+\s*$")
TABLES = {"Table S1.": "S1", "Table S2.": "S2"}
EXCLUDE_TRUTH = [  # (file, cls, feature_types, aliases) already in truth, non-seed
    ("benchmark/truth/dmel/dme_mirgenedb_v3.0.gff3", "miRNA", {"miRNA_primary_transcript", "pre_miRNA"}, None),
    ("benchmark/truth/dmel/dme_mirbase_v23.1.gff3", "miRNA", {"miRNA_primary_transcript", "pre_miRNA"},
     "benchmark/truth/dmel/dme_contig_aliases.tsv"),
    ("benchmark/truth/dmel/dme_piRNAclusterDB2_proTRAC.gff3", "piRNA", None, None),
]
SEED = "benchmark/truth/seed/piC_ref_dm6.bed"


def parse_tables(txt: str) -> list[tuple[str, int, int, str]]:
    rows, table = [], None
    for line in txt.splitlines():
        s = line.strip()
        for key, tag in TABLES.items():
            if s.startswith(key):
                table = tag
        m = ROW.match(line)
        if table and m:
            chrom, start, stop, ann, size = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4).strip(), int(m.group(5))
            assert stop - start == size, line
            rows.append((chrom, start, stop, f"Czech2008_{table}:{ann.replace(' ', '_')}"))
    return rows


def union(rows):
    out = []
    for chrom, s, e, n in sorted(rows):
        if out and out[-1][0] == chrom and s < out[-1][2]:
            c, s0, e0, n0 = out[-1]
            out[-1] = (c, s0, max(e0, e), n0 if n in n0 else f"{n0}+{n.split(':')[0]}:{n.split(':', 1)[1]}")
        else:
            out.append((chrom, s, e, n))
    return out


def overlaps(regs, chrom, s, e):
    return any(r.chrom == chrom and r.start < e and s < r.end for r in regs)


def main(pdf: Path, chain: Path, out: Path, work: Path) -> None:
    txt_path = work / "suppl.txt"
    subprocess.run(["pdftotext", "-layout", str(pdf), str(txt_path)], check=True)
    rows = parse_tables(txt_path.read_text())
    n_parsed = len(rows)
    rows = [r for r in rows if "vector" not in r[3]]
    n_novec = len(rows)
    merged = union(rows)
    src = work / "dm3.bed"
    src.write_text("".join(f"chr{c}\t{s}\t{e}\t{n}\n" for c, s, e, n in merged))
    lifted, unmapped = work / "dm6.bed", work / "dm6.bed.unmap"
    subprocess.run(["CrossMap", "bed", str(chain), str(src), str(lifted)], check=True, stdout=subprocess.DEVNULL)
    ok = [l.split("\t") for l in lifted.read_text().splitlines() if l]
    ok = [(f[0][3:] if f[0].startswith("chr") else f[0], int(f[1]), int(f[2]), f[3]) for f in ok]
    ok = [r for r in ok if r[2] - r[1] >= 2]
    ex = []
    for path, cls, types, al in EXCLUDE_TRUTH:
        ex.append((path, load_regions(REPO / path, cls, path, feature_types=types,
                                      aliases=load_aliases(REPO / al) if al else None)))
    seed = load_regions(REPO / SEED, "piRNA", SEED)
    kept, dropped = [], {p: 0 for p, _ in ex}
    for r in ok:
        hit = next((p for p, regs in ex if overlaps(regs, r[0], r[1], r[2])), None)
        if hit:
            dropped[hit] += 1
        else:
            kept.append(r)
    kept.sort(key=lambda r: (r[0], r[1]))
    out.write_text("".join(f"{c}\t{s}\t{e}\t{n}\t0\t.\n" for c, s, e, n in kept))
    print(f"parsed rows {n_parsed}; without transgene 'vector' rows {n_novec}; union of S1/S2 {len(merged)}; "
          f"lifted {len(ok) + 0} (unmapped/too short {len(merged) - len(ok)}); dropped by overlap {dropped}; "
          f"kept {len(kept)}; kept overlapping tuning seed piC {sum(overlaps(seed, *r[:3]) for r in kept)}", file=sys.stderr)


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    w = Path(sys.argv[4]) if len(sys.argv) > 4 else Path(tempfile.mkdtemp())
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), w)
