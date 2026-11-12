"""`report` subcommand: build the self-contained HTML report from classify output."""
import math
import os
import sys

import pysam

import report_data
import report_html
import report_signals
import report_text
from utils import log


_LARGE_DETAIL = 5000  # above this, warn about file size


def _fragments_path(args):
    """--fragments, else the .fragments.tsv beside --clusters when it exists."""
    if args.fragments:
        return args.fragments
    if args.clusters.endswith(".clusters.tsv"):
        cand = args.clusters[:-len(".clusters.tsv")] + ".fragments.tsv"
        if os.path.isfile(cand):
            return cand
    return None


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
    _need(args.clusters, "clusters TSV")
    _need(args.bam, "BAM")
    _need_index(args.bam)
    frags = _fragments_path(args)
    if frags:
        _need(frags, "fragments TSV")

    clusters = report_data.load_clusters(args.clusters)
    sel = report_data.select(clusters, args.min_reads, args.max_detail)
    log(f"report: {len(sel.shown)} clusters shown, {sel.n_hidden} hidden "
        f"(below {args.min_reads:g} reads) of {sel.n_total}", "INFO")
    if sel.capped:
        log(f"report: detail capped at the top {args.max_detail} of {len(sel.shown)} "
            f"clusters (raise --max-detail to include more)", "WARN")

    if args.max_detail > _LARGE_DETAIL:
        log(f"report: --max-detail {args.max_detail} is above {_LARGE_DETAIL}; each detail "
            f"page is about 7 KB (about 7 MB per 1000 pages), so the file can get large "
            f"and slow to open in a browser", "WARN")

    summary = report_data.library_summary(clusters, frags)
    signals = report_signals.signals_for(args.bam, sel.detail, threads=args.threads)
    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"{args.prefix}.report.html")
    title = args.title or report_text.STRINGS["ui.title.default"]
    report_html.build_report(sel, summary, signals, title, out)
    log(f"report written: {out}", "INFO")
