"""Run RepeatMasker on the classified loci only, not the whole genome.

The output `<prefix>.te.gff3` is in genome coordinates and is meant to be passed
back to classify as `--annotation te:<that file>`.
"""
import csv
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pysam

from utils import log, which

MAX_LOCUS = 1_000_000
DEFAULT_NCRNA_CLASSES = "ambiguous,siRNA,phasiRNA,phasiRNA24,miRNA"


def _die(msg):
    log(msg, "ERROR")
    sys.exit(1)


def write_locus_fasta(clusters_tsv, genome, out_fa, map_path, flank=0,
                      scan_fa=None, scan_classes=None, scan_ids=None):
    """One record per locus, named L1, L2, ... plus a mapping file.

    RepeatMasker limits sequence ids to 50 characters, so `chrom:start-end`
    names do not survive long scaffold names. map_path gets one tab-separated
    row per record: id, chrom, start0 (0-based, of the fetched span), end.
    With flank > 0 the span is padded and the map holds the span fetched.
    out_fa None: no loci FASTA (ncRNA-only mode). scan_fa: a second FASTA with
    the same ids, only for loci whose `cls` is in scan_classes (None = all);
    scan_ids gets {id: cls} for them.
    """
    n = 0
    with pysam.FastaFile(genome) as fa, open(clusters_tsv) as fh, \
            open(out_fa or os.devnull, "w") as out, open(map_path, "w") as mp, \
            open(scan_fa or os.devnull, "w") as scan:
        rd = csv.DictReader(fh, delimiter="\t")
        need = ["chrom", "start", "end"] + (["cls"] if scan_classes is not None else [])
        for col in need:
            if col not in (rd.fieldnames or []):
                _die(f"clusters TSV {clusters_tsv} has no '{col}' column")
        lengths = dict(zip(fa.references, fa.lengths))
        for row in rd:
            chrom, s, e = row["chrom"], int(row["start"]), int(row["end"])
            if chrom not in lengths:
                log(f"locus on unknown contig {chrom}; skipped", "WARN")
                continue
            if e - s > MAX_LOCUS:
                log(f"locus {chrom}:{s}-{e} is {e - s} bp (> 1 Mb); skipped", "WARN")
                continue
            s, e = max(0, s - flank), min(lengths[chrom], e + flank)
            if e <= s:
                continue
            n += 1
            seq = fa.fetch(chrom, s, e)
            out.write(f">L{n}\n{seq}\n")
            mp.write(f"L{n}\t{chrom}\t{s}\t{e}\n")
            if scan_fa and (scan_classes is None or row.get("cls") in scan_classes):
                scan.write(f">L{n}\n{seq}\n")
                if scan_ids is not None:
                    scan_ids[f"L{n}"] = row.get("cls", "")
    return n


def read_loci_map(map_path):
    """{id: (chrom, start0)} from the mapping file."""
    out = {}
    with open(map_path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 3:
                out[cols[0]] = (cols[1], int(cols[2]))
    return out


def rm_out_to_gff(out_path, map_path, gff_path):
    """RepeatMasker .out with L<n> queries -> GFF3 in genome coordinates."""
    loci = read_loci_map(map_path)
    n = 0
    with open(out_path) as fh, open(gff_path, "w") as out:
        out.write("##gff-version 3\n")
        for line in fh:
            cols = line.split()
            if len(cols) < 11:
                continue
            try:
                begin, end = int(cols[5]), int(cols[6])
            except ValueError:
                continue
            if cols[4] not in loci:
                log(f"query name {cols[4]} is not in the loci map; skipped", "WARN")
                continue
            chrom, start0 = loci[cols[4]]
            if begin < 1 or begin > end:
                continue
            strand = "-" if cols[8] == "C" else "+"
            out.write(f"{chrom}\tRepeatMasker\trepeat_region\t{start0 + begin}\t"
                      f"{start0 + end}\t{cols[0]}\t{strand}\t.\t"
                      f"Name={cols[9]};Class={cols[10]}\n")
            n += 1
    return n


def _run_repeatmasker(cmd):
    """The one subprocess seam; tests replace it."""
    subprocess.run(cmd, check=True)


def ncrna_class(family):
    """Coarse RNA type from an Rfam family name."""
    f = family.lower()
    if f.startswith("mir"):
        return "miRNA"
    if f.startswith("sno"):
        return "snoRNA"
    if "trna" in f:
        return "tRNA"
    if any(k in f for k in ("rrna", "lsu", "ssu", "5s", "5_8s")):
        return "rRNA"
    if "snrna" in f or re.match(r"^u\d", f):
        return "snRNA"
    return "ncRNA"


def parse_tblout(path, scan_ids, loci):
    """Hits of a cmscan --fmt 2 tblout with inc '!', as GFF3 rows in genome coordinates.

    Columns (1-based): 2 family, 3 accession, 5 query, 10/11 seq from/to,
    12 strand, 18 E-value, 19 inc. Returns (rows, n_loci_hit).
    """
    rows, hit = [], set()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            if len(c) < 19 or c[18] != "!":
                continue
            q = c[4]
            if q not in loci:
                log(f"query name {q} is not in the loci map; skipped", "WARN")
                continue
            try:
                a, b = sorted((int(c[9]), int(c[10])))
            except ValueError:
                continue
            chrom, start0 = loci[q]
            strand = c[11] if c[11] in ("+", "-") else "."
            rows.append((chrom, start0 + a, start0 + b, strand, c[1], c[2],
                         scan_ids.get(q, ""), c[17]))
            hit.add(q)
    rows.sort(key=lambda r: (r[0], r[1], r[2]))
    return rows, len(hit)


def write_ncrna_gff(rows, gff_path):
    with open(gff_path, "w") as out:
        out.write("##gff-version 3\n")
        for chrom, s, e, strand, fam, acc, cls, ev in rows:
            out.write(f"{chrom}\tInfernal\tncRNA\t{s}\t{e}\t.\t{strand}\t.\t"
                      f"Name={fam};Accession={acc};Class={ncrna_class(fam)};"
                      f"evalue={ev}\n")


def _run_tool(cmd):
    """Subprocess seam for cmpress/cmscan."""
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)


def _z_value(scan_fa):
    bp = 0
    with open(scan_fa) as fh:
        for line in fh:
            if not line.startswith(">"):
                bp += len(line.strip())
    z = max(1.0, 2 * bp / 1e6)
    return f"{z:.3f}".rstrip("0").rstrip(".")


def _parse_classes(value):
    value = (value or DEFAULT_NCRNA_CLASSES).strip()
    if value == "all":
        return None
    return {c.strip() for c in value.split(",") if c.strip()}


def run_annotate_loci(args):
    cm = getattr(args, "ncrna_cm", None)
    te = bool(args.species or args.te_library)
    if args.species and args.te_library:
        _die("give at most one of --species / --te-library" if cm
             else "give exactly one of --species or --te-library")
    if not te and not cm:
        _die("give exactly one of --species or --te-library" if not cm else
             "give --species or --te-library, or --ncrna-cm")
    if te and not which("RepeatMasker"):
        _die("RepeatMasker is not on PATH. Install it with "
             "`conda install -c bioconda repeatmasker` (it is optional for "
             "the rest of srnaprofiler)")
    hint = "Install it with `conda install -c bioconda infernal`"
    if cm:
        if not which("cmscan"):
            _die(f"cmscan is not on PATH. {hint} (it is optional for the rest of srnaprofiler)")
        if not Path(str(cm) + ".i1m").exists() and not which("cmpress"):
            _die(f"cmpress is not on PATH. {hint} (it is optional for the rest of srnaprofiler)")
        if not os.path.isfile(cm):
            _die(f"--ncrna-cm {cm}: no such file")
        clanin = getattr(args, "ncrna_clanin", None)
        if clanin and not os.path.isfile(clanin):
            _die(f"--ncrna-clanin {clanin}: no such file")
    outdir = Path(args.outdir)
    prefix = args.prefix
    fa = outdir / f"{prefix}.loci.fa"
    mp = outdir / f"{prefix}.loci.map.tsv"
    rmdir = outdir / f"{prefix}.repeatmasker"
    gff = outdir / f"{prefix}.te.gff3"
    scan = outdir / f"{prefix}.ncrna.scan.fa"
    tbl = outdir / f"{prefix}.ncrna.tblout"
    ngff = outdir / f"{prefix}.ncrna.gff3"
    for out_path, on in ((gff, te), (ngff, bool(cm))):
        if on and out_path.exists():
            _die(f"{out_path} exists; refusing to overwrite it. Remove it or pick "
                 f"another --prefix/--outdir")
    outdir.mkdir(parents=True, exist_ok=True)
    # Scratch this run creates is removed on failure so a rerun is possible;
    # anything already there before the run is left alone.
    cand = [mp] + ([fa, rmdir] if te else []) + ([scan, tbl] if cm else [])
    created = [p for p in cand if not p.exists()]

    def fail(msg):
        for p in created:
            if p.is_dir():
                shutil.rmtree(p, ignore_errors=True)
            elif p.exists():
                p.unlink()
        _die(msg)

    scan_ids = {}
    classes = _parse_classes(getattr(args, "ncrna_classes", None)) if cm else None
    n = write_locus_fasta(args.clusters, args.genome, str(fa) if te else None, str(mp),
                          scan_fa=str(scan) if cm else None,
                          scan_classes=classes if cm else None, scan_ids=scan_ids)
    if n == 0:
        fail("no loci to annotate")
    if te:
        log(f"{n} loci written to {fa} (ids mapped back by {mp})")
        rmdir.mkdir(exist_ok=True)
        cmd = ["RepeatMasker", "-pa", str(args.threads)]
        cmd += ["-species", args.species] if args.species else ["-lib", args.te_library]
        cmd += ["-dir", str(rmdir), str(fa)]
        try:
            _run_repeatmasker(cmd)
        except (subprocess.CalledProcessError, OSError) as exc:
            fail(f"RepeatMasker failed: {exc}")
        out = rmdir / (fa.name + ".out")
        if not out.exists():
            fail(f"RepeatMasker produced no {out.name} in {rmdir} "
                 f"(no repeats found, or it failed)")
        k = rm_out_to_gff(str(out), str(mp), str(gff))
        log(f"{k} repeat features -> {gff}; use it as --annotation te:{gff}")
    if cm:
        if not scan_ids:
            fail(f"no locus of the classes {sorted(classes)} to scan for structural RNA "
                 f"(use --ncrna-classes all to scan every locus)")
        t0 = time.time()
        try:
            if not Path(str(cm) + ".i1m").exists():
                log(f"{cm}.i1m missing: running cmpress -F {cm}")
                _run_tool(["cmpress", "-F", str(cm)])
            cmd = ["cmscan", "-Z", _z_value(str(scan)), "--cut_ga", "--rfam",
                   "--nohmmonly", "--fmt", "2", "--cpu", str(args.threads)]
            if clanin:
                cmd += ["--clanin", str(clanin), "--oclan", "--oskip"]
            cmd += ["--tblout", str(tbl), "-o", "/dev/null", str(cm), str(scan)]
            log(f"cmscan on {len(scan_ids)} of {n} loci ({'all classes' if classes is None else ','.join(sorted(classes))})")
            _run_tool(cmd)
        except (subprocess.CalledProcessError, OSError) as exc:
            fail(f"cmscan/cmpress failed: {exc}")
        if not tbl.exists():
            fail(f"cmscan produced no {tbl.name}")
        rows, nl = parse_tblout(str(tbl), scan_ids, read_loci_map(str(mp)))
        write_ncrna_gff(rows, str(ngff))
        for p in (scan, tbl):
            p.unlink(missing_ok=True)
        log(f"{len(rows)} structural RNA hits on {nl} loci in {time.time() - t0:.0f} s "
            f"-> {ngff}; use it as --annotation ncrna:{ngff}")
