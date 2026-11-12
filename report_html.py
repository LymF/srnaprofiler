"""Assemble the self-contained HTML report: inline CSS, JS and one JSON data block.

Payload (JSON in <script id="report-data">), read by the front end:
  library, counts, warnings, strings (STRINGS table), rows, details, meta.
  rows[i]: n name, h chrom, s start, e end, r reads, c class, f confidence (null when
    ambiguous), z dominant length, v dot levels (one digit per dot, see _LEVEL),
    d 1 if details[str(i)] exists, nc nearest class and u reason (ambiguous only).
  details[str(i)]:
    v verdict key in strings (verdict.<class> or verdict.generic; the ambiguous
      one takes its reason from the row's u), seq major sequence, rpm,
    strip [[key, level, value|null]]: keys size (value = dominant length), 1u,
      10a (animal) or 5a (plant), strand (value = plus-strand share), pingpong
      (z, animal), phasing (p), multimap (share); level is a _LEVEL digit string,
    ctx [reads_pct, conf_pct|null] (percentiles among shown clusters),
    sig: {"err":1} | {"none":1} | quantised figure data. Every array is a pair
      [max, ints]: ints 0..99, value = int * max / 99, trailing zeros dropped
      (pad with 0 to the axis length), [0, []] when all zero.
      dom dominant length; sz {lo, A|C|G|U: [plus, minus]} reads per length from lo;
      ph {+,-} phasing, offsets 1..100; pp ping-pong, offsets 1..30;
      cv {s start, b bin size, +, -} mean depth per bin; en {b, +, -} 5' ends per bin;
      lg {+,-}: per position [A,C,G,U] ints 0..99 (fraction * 99);
      sq rows [sequence, length, strand, reads, 5' base] exact.
  meta: title, floor, n_total, kingdom, mode, version, github, tagline, install
    (two lines), defining {class: [check keys the class is decided by]}, conf_min,
    conf_k, dot_names, dot_legend, kt, n_hidden, capped, n_detail, reads_shown, summary.
"""
import html
import json
import math
import os
import re

import report_data
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


def _quant(values):
    """[max, ints]: ints 0..99 = round(99 * v / max), trailing zeros dropped.

    value ~ int * max / 99, error at most max / 198. All-zero input is [0, []].
    """
    vals = [v if (v is not None and math.isfinite(v) and v > 0) else 0.0 for v in values]
    m = max(vals, default=0.0)
    if m <= 0:
        return [0, []]
    q = [round(99 * v / m) for v in vals]
    while q and q[-1] == 0:
        q.pop()
    return [_n(float(f"{m:.6g}"), 6), q]


def dequantise(pair, n=None):
    """Inverse of _quant, padded with zeros to n values (test and table-view use)."""
    m, q = pair
    out = [i * m / 99 for i in q]
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
    return {
        "dom": sig.get("dominant"),
        "sz": dict({"lo": sz.get("lo", 18)},
                   **{b: [_quant(sz[b]["+"]), _quant(sz[b]["-"])] for b in "ACGU"}),
        "ph": {s: _quant(sig["phasing"][s]) for s in "+-"},
        "pp": _quant(sig["pingpong"]),
        "cv": {"s": _n(cov["start"], 0), "b": cov["bin"],
               "+": _quant(cov["+"]), "-": _quant(cov["-"])},
        "en": {"b": ends["bin"], "+": _quant(ends["+"]), "-": _quant(ends["-"])},
        "lg": {s: [[round(99 * (p.get(b) or 0)) for b in "ACGU"] for p in sig["logo"][s]]
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
    if amb:
        r["nc"] = rt.nearest_class(c)
    if ex["unresolved"]:
        r["u"] = ex["unresolved"]
    return r


def _detail(c, ex, sig, ctx):
    return {"v": ex["verdict"],
            "strip": [[k, _LEVEL[lv], _n(v, 4)] for k, lv, v in ex["strip"]],
            "ctx": [ctx["reads_pct"], ctx["conf_pct"]],
            "seq": c.major_rna, "rpm": _n(c.reads_rpm), "sig": _signals(sig)}


def _payload(sel, summary, signals, title):
    summary = summary or {}
    detail_ids = {id(c) for c in sel.detail}
    ctx = report_data.library_context(sel.shown)
    rows, details, counts = [], {}, {}
    kingdom = sel.kingdom
    for i, c in enumerate(sel.shown):
        ex = rt.explain(c, kingdom)
        has = id(c) in detail_ids
        rows.append(_row(i, c, has, ex, kingdom))
        counts[c.cls] = counts.get(c.cls, 0) + 1
        if has:
            details[str(i)] = _detail(c, ex, signals.get(c.name), ctx[c.name])
    lib = {k: _n(summary.get(k), 4) for k in
           ("library_reads", "cluster_reads", "unique_share", "fragments", "fragment_reads")}
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
                 "n_detail": len(sel.detail),
                 "reads_shown": _n(sum(c.reads or 0 for c in sel.shown), 1),
                 "summary": rt.summary_sentence(counts, sel.floor)},
    }


def _footer():
    return html.escape(rt.STRINGS["ui.footer"].format(
        present=f"{rt.DISPLAY_PRESENT:g}", weak=f"{rt.DISPLAY_WEAK:g}",
        p10=f"{rt.DISPLAY_10A_PRESENT:g}", w10=f"{rt.DISPLAY_10A_WEAK:g}",
        psize=f"{rt.DISPLAY_SIZE_PRESENT:g}", wsize=f"{rt.DISPLAY_SIZE_WEAK:g}",
        mmod=f"{rt.DISPLAY_MULTIMAP_MODERATE:g}", mhigh=f"{rt.DISPLAY_MULTIMAP_HIGH:g}"))


def build_report(sel, summary, signals, title, out_path):
    data = json.dumps(_payload(sel, summary, signals or {}, title),
                      separators=(",", ":"), allow_nan=False)
    # "<" cannot close the script element; JSON still decodes it back
    data = data.replace("<", "\\u003c")
    parts = {"CSS": _read("report.css"), "JS": _read("app.js"), "DATA": data,
             "TITLE": html.escape(title or rt.STRINGS["ui.title.default"]), "FOOTER": _footer()}
    # one pass: substituted text is never rescanned for tokens
    out = re.sub(r"%%(CSS|JS|DATA|TITLE|FOOTER)%%", lambda m: parts[m.group(1)],
                 _read("template.html"))
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(out)
