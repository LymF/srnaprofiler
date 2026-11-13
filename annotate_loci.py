"""Run RepeatMasker on the classified loci only, not the whole genome.

The output `<prefix>.te.gff3` is in genome coordinates and is meant to be passed
back to classify as `--annotation te:<that file>`.
"""
import csv
import subprocess
import sys
from pathlib import Path

import pysam

from utils import log, which

MAX_LOCUS = 1_000_000


def _die(msg):
    log(msg, "ERROR")
    sys.exit(1)


def write_locus_fasta(clusters_tsv, genome, out_fa, flank=0):
    """One record per locus, named chrom:start-end (0-based start, as in the TSV).

    With flank > 0 the record is padded and named by the span actually
    fetched, so the offsets in rm_out_to_gff stay exact.
    """
    n = 0
    with pysam.FastaFile(genome) as fa, open(clusters_tsv) as fh, \
            open(out_fa, "w") as out:
        rd = csv.DictReader(fh, delimiter="\t")
        for col in ("chrom", "start", "end"):
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
            out.write(f">{chrom}:{s}-{e}\n{fa.fetch(chrom, s, e)}\n")
            n += 1
    return n


def rm_out_to_gff(out_path, gff_path):
    """RepeatMasker .out with chrom:start-end queries -> GFF3 in genome coordinates."""
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
            chrom, _, span = cols[4].rpartition(":")
            try:
                start0 = int(span.split("-")[0])
            except ValueError:
                log(f"query name {cols[4]} is not chrom:start-end; skipped", "WARN")
                continue
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


def run_annotate_loci(args):
    if bool(args.species) == bool(args.te_library):
        _die("give exactly one of --species or --te-library")
    if not which("RepeatMasker"):
        _die("RepeatMasker is not on PATH. Install it with "
             "`conda install -c bioconda repeatmasker` (it is optional for "
             "the rest of srnaprofiler)")
    outdir = Path(args.outdir)
    prefix = args.prefix
    fa = outdir / f"{prefix}.loci.fa"
    rmdir = outdir / f"{prefix}.repeatmasker"
    gff = outdir / f"{prefix}.te.gff3"
    for p in (fa, rmdir, gff):
        if p.exists():
            _die(f"{p} exists; refusing to overwrite it. Remove it or pick "
                 f"another --prefix/--outdir")
    outdir.mkdir(parents=True, exist_ok=True)
    n = write_locus_fasta(args.clusters, args.genome, str(fa))
    if n == 0:
        _die("no loci to annotate")
    log(f"{n} loci written to {fa}")
    rmdir.mkdir()
    cmd = ["RepeatMasker", "-pa", str(args.threads)]
    cmd += ["-species", args.species] if args.species else ["-lib", args.te_library]
    cmd += ["-dir", str(rmdir), str(fa)]
    try:
        _run_repeatmasker(cmd)
    except (subprocess.CalledProcessError, OSError) as exc:
        _die(f"RepeatMasker failed: {exc}")
    out = rmdir / (fa.name + ".out")
    if not out.exists():
        _die(f"RepeatMasker produced no {out.name} in {rmdir} "
             f"(no repeats found, or it failed)")
    k = rm_out_to_gff(str(out), str(gff))
    log(f"{k} repeat features -> {gff}; use it as --annotation te:{gff}")
