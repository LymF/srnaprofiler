"""`report` subcommand: build the self-contained HTML report from classify output."""
import math
import os
import sys

import pysam

import annotation
import report_data
import report_genome
import report_html
import report_mirna
import report_qc
import report_signals
import report_text
import samio
from utils import log


_LARGE_DETAIL = 5000  # above this, warn about file size
_GM_WARN_BYTES = 20 * 1024 * 1024   # warn when the genome view adds more than this
_GM_BINS = 200                      # coverage bins per window in the genome view


def _load_layers(args):
    """Annotation layers for the genome view, or None. A bad label, a missing file or an
    unreadable (truncated/fake .gz) file stops here, before any BAM pass."""
    specs = getattr(args, "annotation", None) or []
    if not specs:
        return None
    try:
        annotation.check_layer_specs(specs)
        return annotation.load_layers(specs, set(samio.bam_references(args.bam)))
    except (ValueError, FileNotFoundError) as e:
        _fail(str(e))


def _genome(args, sel, layers):
    """{cluster name: window dict + "cv"} for the detailed loci, or None without layers."""
    if layers is None:
        return None
    lengths = samio.reference_lengths(args.bam)
    gm = report_genome.build_windows(sel.detail, layers, args.flank, lengths)
    wins = {name: (c.chrom, g["w"][0], g["w"][1])
            for c in sel.detail for name, g in [(c.name, gm[c.name])]}
    cv = report_signals.window_tracks_for(args.bam, wins, threads=args.threads,
                                          max_bins=_GM_BINS)
    for name, g in gm.items():
        g["cv"] = cv[name]
    return gm


def _fragments_path(args):
    """--fragments, else the .fragments.tsv beside --clusters when it exists."""
    if args.fragments:
        return args.fragments
    if args.clusters.endswith(".clusters.tsv"):
        cand = args.clusters[:-len(".clusters.tsv")] + ".fragments.tsv"
        if os.path.isfile(cand):
            return cand
    return None


def _mirna_path(args):
    """--mirna, else the .mirna.tsv beside --clusters when it exists."""
    if getattr(args, "mirna", None):
        return args.mirna
    if args.clusters.endswith(".clusters.tsv"):
        cand = args.clusters[:-len(".clusters.tsv")] + ".mirna.tsv"
        if os.path.isfile(cand):
            return cand
    return None


def _mirna_figures(path, sel):
    """({cluster name: (Row, fold)}, {cluster name: hairpin spec}) for the shown miRNA clusters
    that have a detail page and a row in the miRNA table."""
    try:
        rows = report_mirna.load(path)
    except ValueError as e:
        _fail(str(e))
    hit = report_mirna.match(sel.detail, rows)
    missing = report_mirna.fold_tool_missing()
    mirna, specs = {}, {}
    for name, row in hit.items():
        db, mfe = report_mirna.fold(row.precursor_seq)
        mirna[name] = (row, (db, mfe, missing, report_mirna.layout(db)))
        spec = report_mirna.hairpin_spec(row)
        if spec:
            specs[name] = spec
    n_cl = sum(1 for c in sel.detail if c.cls == "miRNA")
    log(f"report: {len(hit)} of {n_cl} miRNA clusters matched to {os.path.basename(path)} "
        f"({len(rows)} rows)", "INFO")
    if hit and missing and any(m[1][0] is None for m in mirna.values()):
        log("report: RNAfold not found (executable or python module RNA): the hairpin "
            "figures omit the structure", "WARN")
    return mirna, specs


def _need(path, what):
    if not os.path.isfile(path):
        log(f"{what} not found: {path}", "ERROR")
        sys.exit(f"report: {what} not found: {path}")


def _fail(msg):
    log(msg, "ERROR")
    sys.exit(f"report: {msg}")


def _need_index(bam):
    """Every figure fetches a region: without an index each fetch would fail."""
    try:
        with pysam.AlignmentFile(bam, "rb") as fh:
            ok = fh.check_index()
    except (ValueError, OSError) as e:
        _fail(f"BAM index not usable for {bam} ({e}); create it with `samtools index`")
    if not ok:
        _fail(f"BAM index not found for {bam}; create it with `samtools index`")


def run_report(args):
    if args.max_detail < 0:
        _fail(f"--max-detail must be 0 or more (got {args.max_detail})")
    if not math.isfinite(args.min_reads) or args.min_reads < 0:
        _fail(f"--min-reads must be a finite number, 0 or more (got {args.min_reads})")
    if args.threads < 1:
        _fail(f"--threads must be 1 or more (got {args.threads})")
    if getattr(args, "flank", 2000) < 0:
        _fail(f"--flank must be 0 or more (got {args.flank})")
    _need(args.clusters, "clusters TSV")
    _need(args.bam, "BAM")
    _need_index(args.bam)
    layers = _load_layers(args)   # validated before any BAM pass
    frags = _fragments_path(args)
    if frags:
        _need(frags, "fragments TSV")
    mirna_path = _mirna_path(args)
    if mirna_path:
        _need(mirna_path, "miRNA TSV")
    if getattr(args, "fastp", None):
        _need(args.fastp, "fastp JSON")
    if getattr(args, "qc_summary", None):
        _need(args.qc_summary, "QC summary")

    clusters = report_data.load_clusters(args.clusters)
    sel = report_data.select(clusters, args.min_reads, args.max_detail)
    log(f"report: {len(sel.shown)} clusters shown, {sel.n_hidden} hidden "
        f"(below {args.min_reads:g} reads) of {sel.n_total}", "INFO")
    if sel.capped:
        log(f"report: detail capped at {args.max_detail} of {len(sel.shown)} "
            f"clusters, each class keeping a share (raise --max-detail to include more)", "WARN")

    if args.max_detail > _LARGE_DETAIL:
        log(f"report: --max-detail {args.max_detail} is above {_LARGE_DETAIL}; each detail "
            f"page is about 7 KB (about 7 MB per 1000 pages), so the file can get large "
            f"and slow to open in a browser", "WARN")

    summary = report_data.library_summary(clusters, frags)
    if not getattr(args, "no_library_qc", False):
        lengths = report_signals.library_lengths(args.bam, threads=args.threads)
        fastp, qcs = report_qc.trim_paths(args.clusters, getattr(args, "fastp", None),
                                          getattr(args, "qc_summary", None))
        summary["qc"] = report_qc.build(lengths, report_qc.load_trim(fastp, qcs))
    mirna, specs = _mirna_figures(mirna_path, sel) if mirna_path else ({}, {})
    signals = report_signals.signals_for(args.bam, sel.detail, threads=args.threads,
                                         hairpins=specs)
    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"{args.prefix}.report.html")
    title = args.title or report_text.STRINGS["ui.title.default"]
    genome = _genome(args, sel, layers)
    report_html.build_report(sel, summary, signals, title, out, mirna=mirna, genome=genome,
                             flank=getattr(args, "flank", None))
    if genome:
        added = sum(len(report_html.gm_json(g)) for g in genome.values())
        log(f"report: genome view for {len(genome)} windows, adds {added / 1e6:.2f} MB "
            f"({added // len(genome)} bytes per locus)", "INFO")
        if added > _GM_WARN_BYTES:
            log(f"report: the genome view adds {added / 1e6:.1f} MB (over 20 MB); lower "
                f"--max-detail or --flank to shrink the file", "WARN")
    log(f"report written: {out}", "INFO")
