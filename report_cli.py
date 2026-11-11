"""`report` subcommand: build the self-contained HTML report from classify output."""
import os
import sys

import report_data
import report_html
import report_signals
import report_text
from utils import log


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


def run_report(args):
    _need(args.clusters, "clusters TSV")
    _need(args.bam, "BAM")
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

    summary = report_data.library_summary(clusters, frags)
    signals = report_signals.signals_for(args.bam, sel.detail)
    os.makedirs(args.outdir, exist_ok=True)
    out = os.path.join(args.outdir, f"{args.prefix}.report.html")
    title = args.title or report_text.STRINGS["ui.title.default"]
    report_html.build_report(sel, summary, signals, title, out)
    log(f"report written: {out}", "INFO")
