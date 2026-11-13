"""Run RepeatMasker on the classified loci only, not the whole genome.

The output `<prefix>.te.gff3` is in genome coordinates and is meant to be passed
back to classify as `--annotation te:<that file>`.
"""
import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pysam

from utils import log, which

MAX_LOCUS = 1_000_000


def _die(msg):
    log(msg, "ERROR")
    sys.exit(1)


def write_locus_fasta(clusters_tsv, genome, out_fa, map_path, flank=0):
    """One record per locus, named L1, L2, ... plus a mapping file.

    RepeatMasker limits sequence ids to 50 characters, so `chrom:start-end`
    names do not survive long scaffold names. map_path gets one tab-separated
    row per record: id, chrom, start0 (0-based, of the fetched span), end.
    With flank > 0 the span is padded and the map holds the span fetched.
    """
    n = 0
    with pysam.FastaFile(genome) as fa, open(clusters_tsv) as fh, \
            open(out_fa, "w") as out, open(map_path, "w") as mp:
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
            n += 1
            out.write(f">L{n}\n{fa.fetch(chrom, s, e)}\n")
            mp.write(f"L{n}\t{chrom}\t{s}\t{e}\n")
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
    mp = outdir / f"{prefix}.loci.map.tsv"
    rmdir = outdir / f"{prefix}.repeatmasker"
    gff = outdir / f"{prefix}.te.gff3"
    if gff.exists():
        _die(f"{gff} exists; refusing to overwrite it. Remove it or pick "
             f"another --prefix/--outdir")
    outdir.mkdir(parents=True, exist_ok=True)
    # Scratch this run creates is removed on failure so a rerun is possible;
    # anything already there before the run is left alone.
    created = [p for p in (fa, mp, rmdir) if not p.exists()]

    def fail(msg):
        for p in created:
            if p.is_dir():
                shutil.rmtree(p, ignore_errors=True)
            elif p.exists():
                p.unlink()
        _die(msg)

    n = write_locus_fasta(args.clusters, args.genome, str(fa), str(mp))
    if n == 0:
        fail("no loci to annotate")
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
