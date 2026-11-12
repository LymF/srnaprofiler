"""Assemble the self-contained HTML report: inline CSS, JS and one JSON data block.

Payload (JSON in <script id="report-data">), read by the front end:
  library, counts, warnings, strings (STRINGS table), rows, details, meta.
  library.qc (absent with --no-library-qc or without reads): the library quality strip from
    report_qc.build: c 27 weighted read counts (<16, 16..40, >40), u unmapped records,
    fl [[check key, value, status 1 ok / 0 check]], tr {b before, a after, ad adapter-trimmed}.
  rows[i]: n name, h chrom, s start, e end, r reads, c class, f confidence (null when
    ambiguous), z dominant length, v dot levels (one digit per dot, see _LEVEL),
    rv 1 if the classifier flagged the locus for review (else absent),
    d 1 if details[str(i)] exists, nc nearest class and u reason (ambiguous only).
  details[str(i)]:
    v verdict key in strings (verdict.<class> or verdict.generic; the ambiguous
      one takes its reason from the row's u), seq major sequence, rpm,
    strip [[key, level, value|null]]: keys size (value = dominant length), 1u,
      10a (animal) or 5a (plant), strand (value = plus-strand share), pingpong
      (z, animal), phasing (p), multimap (share); level is a _LEVEL digit string,
    ctx [reads_pct, conf_pct|null] (percentiles among shown clusters),
    ck [class, [status, observed text], ...] class checklist, rows in meta.cl[class] order
      (status: 1 met, 0 not met, 2 not available, 3 not tested; class is the nearest class
      for an ambiguous cluster; absent for miRNA, other and ML-ambiguous), see report_checks,
    sig: {"err":1} | {"none":1} | quantised figure data:
      n total aligned reads the panels are drawn from (int), nd weighted reads at
      the dominant length (logo and sizes panel), dom dominant length, nb number
      of bins of cv and en (the arrays below are padded with 0 to their fixed length).
      Pairs [max, ints]: ints 0..99, value = int * max / 99, trailing zeros dropped,
      [0, []] when all zero (error <= max / 198).
      sz {lo, max, A|C|G|U: [plus, minus]}: ONE shared max for the 8 arrays, ints
        0..999, value = int * max / 999 (a stack of 4 bases errs <= 0.2% of max);
        18 slots per array, lengths lo..35.
      ph {+,-} phasing, 100 slots (offsets 1..100); pp ping-pong, 30 slots (1..30);
      cv {s start, b bin size, +, -} mean depth per bin, nb slots;
      en {b, +, -} 5' ends per bin, nb slots;
      lg {+,-}: 26 positions, each [A,C,G,U] ints 0..99 (fraction * 99);
      sq rows [sequence, length, strand, reads, 5' base] exact.
    mi (only for a miRNA cluster matched to a row of --mirna): the hairpin figure and its
      numbers, see report_mirna.figure_data: n name, st strand, L length, seq precursor (RNA),
      m/s/lp [a, b] mature, star and loop (0-based half-open, strand-oriented), db dot-bracket
      and mfe or ns (m RNAfold missing, f failed), rd [[a, b, sequence (lower case matches),
      reads, kind m|s|""]], cov per-base depth, more/mr rows and reads not drawn, anti/ovh
      reads not shown, lo genomic start of the precursor window, nowin 1 when the window could not be placed, t numbers [[key, value]],
      pc plant criteria [[key, status 1/0/2, observed, rule]]. Exact, not quantised.
  meta: title, floor, n_total, kingdom, mode, version, github, tagline, install
    (two lines), defining {class: [check keys the class is decided by]}, conf_min,
    conf_k, dot_names, dot_legend, kt, n_hidden, capped, n_detail, reads_shown, summary,
    cl {class: [[check name key, rule text], ...]} the checklist rows per class.
"""
import base64
import html
import json
import math
import os
import re

import report_checks
import report_data
import report_mirna
import report_text as rt
import version

_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "report_assets")
_LEVEL = {"absent": "0", "weak": "1", "present": "2", "unavailable": "3",
          "low": "4", "moderate": "5", "high": "6"}


def _read(name):
    with open(os.path.join(_DIR, name), encoding="utf-8") as fh:
        return fh.read()


def _n(x, nd=2):
    """Compact number: None for missing/non-finite, int when integral."""
    if x is None:
        return None
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(x):
        return None
    if x == int(x):
        return int(x)
    return round(x, nd)


def _fin(v):
    return v if (v is not None and math.isfinite(v) and v > 0) else 0.0


def _ints(vals, m, top):
    q = [round(top * v / m) for v in vals]
    while q and q[-1] == 0:
        q.pop()
    return q


def _qmax(m):
    # 6 significant digits keep the scale of tiny maxima
    return float(f"{m:.6g}")


def _quant(values):
    """[max, ints]: ints 0..99 = round(99 * v / max), trailing zeros dropped.

    value ~ int * max / 99, error at most max / 198. All-zero input is [0, []].
    """
    vals = [_fin(v) for v in values]
    m = max(vals, default=0.0)
    if m <= 0:
        return [0, []]
    return [_qmax(m), _ints(vals, m, 99)]


def _quant_group(arrays, top=999):
    """(max, [ints...]) with one shared max over all arrays (stacked panels)."""
    arrs = [[_fin(v) for v in a] for a in arrays]
    m = max((max(a, default=0.0) for a in arrs), default=0.0)
    if m <= 0:
        return 0, [[] for _ in arrs]
    return _qmax(m), [_ints(a, m, top) for a in arrs]


def dequantise(pair, n=None, top=99):
    """Inverse of _quant (top=99) or of one sz array with its shared max (top=999)."""
    m, q = pair
    out = [i * m / top for i in q]
    return out + [0.0] * (n - len(out)) if n else out


def _signals(sig):
    """Quantised figure data of one cluster; see the module docstring for the keys."""
    if not sig:
        return None
    # explicit markers instead of empty axes
    if sig.get("failed"):
        return {"err": 1}
    if not sig.get("n_alns"):
        return {"none": 1}
    cov, ends, sz = sig["cov"], sig["ends"], sig["sizes"]
    lo, dom = sz.get("lo", 18), sig.get("dominant")
    m, q = _quant_group([sz[b][s] for b in "ACGU" for s in "+-"])
    nd = 0.0
    if dom is not None and 0 <= dom - lo < len(sz["A"]["+"]):
        nd = sum(_fin(sz[b][s][dom - lo]) for b in "ACGU" for s in "+-")
    return {
        "n": round(sig["n_alns"]), "nd": _n(nd, 1), "dom": dom, "nb": len(cov["+"]),
        "sz": dict({"lo": lo, "max": m},
                   **{b: [q[2 * i], q[2 * i + 1]] for i, b in enumerate("ACGU")}),
        "ph": {s: _quant(sig["phasing"][s]) for s in "+-"},
        "pp": _quant(sig["pingpong"]),
        "cv": {"s": _n(cov["start"], 0), "b": cov["bin"],
               "+": _quant(cov["+"]), "-": _quant(cov["-"])},
        "en": {"b": ends["bin"], "+": _quant(ends["+"]), "-": _quant(ends["-"])},
        "lg": {s: [[round(99 * _fin(p.get(b))) for b in "ACGU"] for p in sig["logo"][s]]
               for s in "+-"},
        "sq": [[q[0], q[1], q[2], _n(q[3], 3), q[4]] for q in sig.get("seqs", [])],
    }


def _row(i, c, has_detail, ex, kingdom):
    st = rt.evidence_status(c, kingdom)
    v = "".join(_LEVEL[x] for x in st.values())
    amb = c.cls == "ambiguous"
    # an ambiguous row's stored confidence is the nearest class's: not shown
    r = {"n": c.name, "h": c.chrom, "s": _n(c.start, 0), "e": _n(c.end, 0), "r": _n(c.reads, 1),
         "c": c.cls, "f": None if amb else _n(c.confidence), "z": c.dicer_call, "v": v,
         "d": 1 if has_detail else 0}
    if c.review:
        r["rv"] = 1
    if amb:
        r["nc"] = rt.nearest_class(c)
    if ex["unresolved"]:
        r["u"] = ex["unresolved"]
    return r


def _detail(c, ex, sig, ctx, kingdom="animal", mirna=None):
    d = {"v": ex["verdict"],
         "strip": [[k, _LEVEL[lv], _n(v, 4)] for k, lv, v in ex["strip"]],
         "ctx": [ctx["reads_pct"], ctx["conf_pct"]],
         "seq": c.major_rna, "rpm": _n(c.reads_rpm), "sig": _signals(sig)}
    ck = report_checks.checklist(c, kingdom)
    if ck:
        d["ck"] = ck
    if mirna:
        row, folded = mirna
        d["mi"] = report_mirna.figure_data(row, folded, (sig or {}).get("hp"))
    return d


def _payload(sel, summary, signals, title, mirna=None):
    summary = summary or {}
    detail_ids = {id(c) for c in sel.detail}
    ctx = report_data.library_context(sel.shown)
    rows, details, counts = [], {}, {}
    kingdom = sel.kingdom
    mirna = mirna or {}
    for i, c in enumerate(sel.shown):
        ex = rt.explain(c, kingdom)
        has = id(c) in detail_ids
        rows.append(_row(i, c, has, ex, kingdom))
        counts[c.cls] = counts.get(c.cls, 0) + 1
        if has:
            details[str(i)] = _detail(c, ex, signals.get(c.name), ctx[c.name], kingdom,
                                      mirna.get(c.name))
    lib = {k: _n(summary.get(k), 4) for k in
           ("library_reads", "cluster_reads", "unique_share", "fragments", "fragment_reads")}
    if summary.get("qc"):
        lib["qc"] = summary["qc"]
    S = rt.STRINGS
    return {
        "library": lib, "counts": counts, "warnings": rt.warnings(sel, summary),
        "rows": rows, "details": details, "strings": S,
        "meta": {"title": title, "floor": _n(sel.floor), "n_total": sel.n_total,
                 "kingdom": kingdom, "mode": sel.mode,
                 "version": version.__version__, "github": S["ui.about.github"],
                 "tagline": S["ui.about.tagline"],
                 "install": [S["ui.about.install.1"], S["ui.about.install.2"]],
                 "defining": {k: sorted(v["defining"])
                              for k, v in report_data.class_specs(kingdom).items()},
                 "conf_min": _n(rt.confidence_floor(kingdom), 3),
                 "conf_k": round(1 / rt.confidence_floor(kingdom)),
                 "dot_names": [S[f"ui.dot.{k}"] for k in rt.DOTS[kingdom]],
                 "dot_legend": rt.dot_legend(kingdom), "kt": rt.kingdom_text(kingdom),
                 "n_hidden": sel.n_hidden, "capped": bool(sel.capped),
                 "n_detail": len(sel.detail), "cl": report_checks.skeleton(kingdom),
                 "reads_shown": _n(sum(c.reads or 0 for c in sel.shown), 1),
                 "summary": rt.summary_sentence(counts, sel.floor)},
    }


# IBM Plex, latin subset (OFL, see report_assets/fonts/LICENSE-IBM-Plex.txt)
_FONTS = (("IBM Plex Sans", "ibm-plex-sans-latin-{w}-normal.woff2"),
          ("IBM Plex Mono", "ibm-plex-mono-latin-{w}-normal.woff2"))


def _font_faces():
    """@font-face rules with the woff2 files inlined as base64 data URIs."""
    out = []
    for family, pattern in _FONTS:
        for w in (400, 500):
            with open(os.path.join(_DIR, "fonts", pattern.format(w=w)), "rb") as fh:
                b64 = base64.b64encode(fh.read()).decode("ascii")
            out.append(f"@font-face{{font-family:'{family}';font-style:normal;font-weight:{w};"
                       f"font-display:swap;src:url(data:font/woff2;base64,{b64}) format('woff2');}}")
    return "\n".join(out) + "\n"


# display-only cutoffs, worded for the checks each kingdom shows
_FOOTER = {"animal": "ui.footer.animal", "plant": "ui.footer.plant"}


def _footer(kingdom="animal"):
    return html.escape(rt.STRINGS[_FOOTER.get(kingdom, "ui.footer.animal")].format(
        present=f"{rt.DISPLAY_PRESENT:g}", weak=f"{rt.DISPLAY_WEAK:g}",
        p10=f"{rt.DISPLAY_10A_PRESENT:g}", w10=f"{rt.DISPLAY_10A_WEAK:g}",
        psize=f"{rt.DISPLAY_SIZE_PRESENT:g}", wsize=f"{rt.DISPLAY_SIZE_WEAK:g}",
        mmod=f"{rt.DISPLAY_MULTIMAP_MODERATE:g}", mhigh=f"{rt.DISPLAY_MULTIMAP_HIGH:g}"))


def build_report(sel, summary, signals, title, out_path, mirna=None):
    """mirna: {cluster name: (report_mirna.Row, (dot-bracket, mfe, tool missing))} or None."""
    data = json.dumps(_payload(sel, summary, signals or {}, title, mirna),
                      separators=(",", ":"), allow_nan=False)
    # "<" cannot close the script element; JSON still decodes it back
    data = data.replace("<", "\\u003c")
    parts = {"CSS": _font_faces() + _read("report.css"),"JS": _read("app.js"), "DATA": data,
             "TITLE": html.escape(title or rt.STRINGS["ui.title.default"]),
             "FOOTER": _footer(sel.kingdom)}
    # one pass: substituted text is never rescanned for tokens
    out = re.sub(r"%%(CSS|JS|DATA|TITLE|FOOTER)%%", lambda m: parts[m.group(1)],
                 _read("template.html"))
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(out)
