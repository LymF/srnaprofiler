"""Per-class plotting: read the classifier TSV + BAM, plot signals grouped by class."""
import os
import shutil
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass

from utils import log

import plots
import class_signals as cs
import class_summary
from samio import read_region, prefers_xm_as_nh

_REQUIRED = ("chrom", "start", "end", "name", "reads", "dicer_call", "cls")


@dataclass
class Locus:
    chrom: str
    start: int
    end: int
    name: str
    cls: str
    reads: float
    dicer_call: str


def load_clusters(tsv_path):
    """Parse the reporter TSV into Locus rows, skipping rows with a blank `cls`."""
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        missing = [c for c in _REQUIRED if c not in idx]
        if missing:
            raise ValueError(f"clusters TSV missing columns: {missing}")
        loci = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cls = f[idx["cls"]]
            if not cls:
                continue
            loci.append(Locus(
                chrom=f[idx["chrom"]], start=int(f[idx["start"]]),
                end=int(f[idx["end"]]), name=f[idx["name"]],
                cls=cls, reads=float(f[idx["reads"]]),
                dicer_call=f[idx["dicer_call"]]))
    return loci


def group_by_class(loci, classes=None):
    """Group Locus rows by `cls`; restrict to `classes` if given.

    Class matching is case-insensitive (the classifier writes "siRNA"/"piRNA",
    but users naturally type "sirna,pirna"); the original `cls` is kept as the
    group key so plot names use the canonical casing.
    """
    wanted = {c.lower() for c in classes} if classes else None
    grouped = defaultdict(list)
    for l in loci:
        if wanted is None or l.cls.lower() in wanted:
            grouped[l.cls].append(l)
    return dict(grouped)


def top_n_loci(loci, n):
    """Return the `n` loci with the largest `reads`, descending; `n <= 0` -> []."""
    if n <= 0:
        return []
    return sorted(loci, key=lambda l: l.reads, reverse=True)[:n]


def dominant_dicer(loci, default=21):
    """Return the most common integer `dicer_call` among loci, ignoring 'N'/''."""
    calls = Counter(int(l.dicer_call) for l in loci if l.dicer_call not in ("N", ""))
    return calls.most_common(1)[0][0] if calls else default


def parse_region(spec, name=None):
    """Parse 'chrom:start-end' (0-based half-open) into a Locus."""
    try:
        chrom, coords = spec.rsplit(":", 1)
        start_s, end_s = coords.split("-")
        start, end = int(start_s), int(end_s)
    except ValueError:
        raise ValueError(f"malformed region (expected chrom:start-end): {spec!r}")
    return Locus(chrom=chrom, start=start, end=end, name=name or spec,
                 cls="", reads=0.0, dicer_call="N")


def run_weblogo(fasta_path, out_pdf, title):
    """Run the external WebLogo 3 CLI; skip (return False) if unavailable/empty."""
    if shutil.which("weblogo") is None:
        log("weblogo CLI not found in PATH; skipping logo", "WARN")
        return False
    if not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        log(f"empty FASTA, skipping logo: {fasta_path}", "WARN")
        return False
    cmd = ["weblogo", "-F", "pdf", "-A", "rna", "-U", "bits", "-s", "large",
           "--resolution", "300", "--errorbars", "NO", "-c", "classic",
           "-P", title, "-f", fasta_path, "-o", out_pdf]
    r = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if r.returncode != 0:
        log(f"weblogo failed: {r.stderr}", "ERROR")
        return False
    return True


def _fetch(bam, locus, xm):
    return list(read_region(bam, locus.chrom, locus.start, locus.end, xm_as_nh=xm))


def _safe_name(name):
    """Filesystem-safe unit name (drvfs/Windows forbids ':')."""
    return name.replace(":", "_").replace("-", "_")


def _write_pingpong_tsv(counts, path, hi=30):
    with open(path, "w") as fh:
        fh.write("overlap\tpairs\n")
        for o in range(1, hi + 1):
            fh.write(f"{o}\t{counts.get(o, 0.0)}\n")


def _write_distance_tsv(counts, path, maxdist=100):
    with open(path, "w") as fh:
        fh.write("distance\tcount\n")
        for d in range(1, maxdist + 1):
            fh.write(f"{d}\t{counts.get(d, 0.0)}\n")


def _write_coverage_tsv(profile, path):
    with open(path, "w") as fh:
        fh.write("bin\tplus\tminus\n")
        for i, (p, m) in enumerate(zip(profile["+"], profile["-"])):
            fh.write(f"{i}\t{p}\t{m}\n")


def _write_pattern_matrix(rows, path):
    with open(path, "w") as fh:
        fh.write("Chromosome#" + "#".join(str(i) for i in range(42)) + "\n")
        for name, row in rows:
            fh.write(name + "#" + "#".join(f"{z:.6f}" for z in row) + "\n")


def _signal_plots_for(name, members, period, size_lo, size_hi, trim,
                      tables, plotsdir, prefix, rpath):
    """Write tables + render all signal plots for one pooled unit (class or locus).

    members: list of (alns, region_start, region_end).
    """
    alns_list = [m[0] for m in members]
    all_alns = [a for alns in alns_list for a in alns]

    if not all_alns:
        log(f"no reads for unit {name}; skipping signal plots", "WARN")
        return

    comp = cs.base_composition(all_alns)
    sd_tsv = tables / f"{prefix}.{name}.sizedist.tsv"
    plots.write_base_composition_tsv(comp, sd_tsv)
    plots.generate_base_distribution_plot(
        sd_tsv, plotsdir / f"{prefix}.{name}.sizedist", rpath, "per_base")

    rows = [(f"{name}_{i}", cs.pattern_row(alns)) for i, alns in enumerate(alns_list)]
    _write_pattern_matrix(rows, tables / f"{prefix}.{name}.pattern.matrix")
    plots.generate_pattern_plot(
        tables / f"{prefix}.{name}.pattern.matrix",
        plotsdir / f"{prefix}.{name}.pattern", rpath)

    pp = Counter()
    for alns in alns_list:
        pp.update(cs.pingpong_offsets(alns, hi=30))
    _write_pingpong_tsv(pp, tables / f"{prefix}.{name}.pingpong.tsv")
    plots.run_r_script("plotPingPong.R",
        [str(tables / f"{prefix}.{name}.pingpong.tsv"), str(plotsdir / f"{prefix}.{name}")],
        rpath)

    for strand, tag in (("+", "plus"), ("-", "minus")):
        ph = Counter()
        for alns in alns_list:
            ph.update(cs.phasing_series(alns, strand, maxdist=100))
        _write_distance_tsv(ph, tables / f"{prefix}.{name}.phasing.{tag}.tsv")
    plots.run_r_script("plotPhasing.R", [
        str(tables / f"{prefix}.{name}.phasing.plus.tsv"),
        str(tables / f"{prefix}.{name}.phasing.minus.tsv"),
        str(period), str(plotsdir / f"{prefix}.{name}")], rpath)

    acc = {"+": [0.0] * 100, "-": [0.0] * 100}
    for alns, rs, re in members:
        prof = cs.coverage_profile(alns, rs, re, bins=100)
        for s in ("+", "-"):
            for i in range(100):
                acc[s][i] += prof[s][i]
    if members:
        for s in ("+", "-"):
            acc[s] = [v / len(members) for v in acc[s]]
    _write_coverage_tsv(acc, tables / f"{prefix}.{name}.coverage.tsv")
    plots.run_r_script("plotCoverageProfile.R",
        [str(tables / f"{prefix}.{name}.coverage.tsv"), str(plotsdir / f"{prefix}.{name}")],
        rpath)

    # logos (sense + antisense) via WebLogo
    for strand, title in (("+", "sense"), ("-", "antisense")):
        seqs = cs.strand_sequences(all_alns, strand, size_lo, size_hi, trim)
        fa = tables / f"{prefix}.{name}.logo.{title}.fasta"
        cs.write_fasta(seqs, fa)
        run_weblogo(str(fa), str(plotsdir / f"{prefix}.{name}.logo.{title}.pdf"), title)


def run_classplots(args):
    from pathlib import Path
    out = Path(args.outdir)
    tables = out / "tables"; plotsdir = out / "plots"
    tables.mkdir(parents=True, exist_ok=True)
    plotsdir.mkdir(parents=True, exist_ok=True)

    if not args.no_summary:
        class_summary.write_summary(args.clusters, str(plotsdir), args.prefix)

    if args.no_signals:
        return

    loci = load_clusters(args.clusters)
    classes = set(args.classes.split(",")) if args.classes else None
    grouped = group_by_class(loci, classes)
    if not grouped:
        available = sorted({l.cls for l in loci})
        log(f"no loci to plot (classes={args.classes!r} matched nothing); "
            f"available classes: {available}", "WARN")
    xm = prefers_xm_as_nh(args.bam)
    rpath = plots.get_r_scripts_path()

    for cls, cls_loci in grouped.items():
        period = dominant_dicer(cls_loci)
        trim = args.logo_trim or period
        members = [(_fetch(args.bam, l, xm), l.start, l.end) for l in cls_loci]
        _signal_plots_for(cls, members, period, args.min_length, args.max_length,
                          trim, tables, plotsdir, args.prefix, rpath)
        for l in top_n_loci(cls_loci, args.top_n):
            m = [(_fetch(args.bam, l, xm), l.start, l.end)]
            _signal_plots_for(f"{cls}.locus.{_safe_name(l.name)}", m, period, args.min_length,
                              args.max_length, trim, tables, plotsdir, args.prefix, rpath)

    # explicit user-requested coordinates (regardless of class)
    default_period = dominant_dicer(loci)
    for spec in (getattr(args, "region", None) or []):
        reg = parse_region(spec)
        trim = args.logo_trim or default_period
        m = [(_fetch(args.bam, reg, xm), reg.start, reg.end)]
        safe = _safe_name(reg.name)
        _signal_plots_for(f"region.{safe}", m, default_period, args.min_length,
                          args.max_length, trim, tables, plotsdir, args.prefix, rpath)
