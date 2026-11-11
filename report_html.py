"""Assemble the self-contained HTML report: inline CSS, JS and one JSON data block."""
import html
import json
import math
import os
import re

import report_text as rt

_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "report_assets")
_LEVEL = {"absent": "0", "weak": "1", "present": "2", "unavailable": "3"}


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


def _arr(xs, nd=2):
    return [_n(v, nd) or 0 for v in xs]


def _signals(sig):
    if not sig:
        return None
    nt = [[round(100 * (p.get(b) or 0)) for b in "ACGU"] for p in sig.get("nt", [])]
    cov = sig.get("coverage", {})
    sz = sig.get("sizes", {})
    return {"sz": [sz.get("lo", 18), _arr(sz.get("plus", [])), _arr(sz.get("minus", []))],
            "nt": nt,
            "cv": [_n(cov.get("start"), 0) or 0, cov.get("bin", 1),
                   _arr(cov.get("plus", []), 1), _arr(cov.get("minus", []), 1)],
            "pp": _arr(sig.get("pingpong", []))}


def _row(i, c, has_detail, ex):
    st = rt.evidence_status(c)
    v = "".join(_LEVEL[st[k]] for k in ("size", "1u", "10a", "pingpong", "phasing"))
    r = {"n": c.name, "h": c.chrom, "s": _n(c.start, 0), "e": _n(c.end, 0), "r": _n(c.reads, 1),
         "c": c.cls, "f": _n(c.confidence), "z": c.dicer_call, "v": v, "d": 1 if has_detail else 0}
    if ex["missing"]:
        r["m"] = ex["missing"]
    return r


def _detail(c, ex, sig):
    return {"verdict": ex["verdict"],
            "checks": [[k["key"], _LEVEL[k["status"]], k["detail"]] for k in ex["checks"]],
            "why": [[w["cls"], w["reason"]] for w in ex["why_not"]],
            "cautions": ex["cautions"], "seq": c.major_rna, "rpm": _n(c.reads_rpm),
            "sig": _signals(sig)}


def _payload(sel, summary, signals, title):
    summary = summary or {}
    detail_ids = {id(c) for c in sel.detail}
    rows, details, counts = [], {}, {}
    for i, c in enumerate(sel.shown):
        ex = rt.explain(c)
        has = id(c) in detail_ids
        rows.append(_row(i, c, has, ex))
        counts[c.cls] = counts.get(c.cls, 0) + 1
        if has:
            details[str(i)] = _detail(c, ex, signals.get(c.name))
    lib = {k: _n(summary.get(k), 4) for k in
           ("library_reads", "cluster_reads", "unique_share", "fragments", "fragment_reads")}
    return {
        "library": lib, "counts": counts, "warnings": rt.warnings(sel, summary),
        "rows": rows, "details": details, "strings": rt.STRINGS, "glossary": rt.GLOSSARY,
        "references": [{"key": k, "text": t} for k, t in rt.REFERENCES.items()],
        "meta": {"title": title, "floor": sel.floor, "n_total": sel.n_total,
                 "n_hidden": sel.n_hidden, "capped": bool(sel.capped),
                 "n_detail": len(sel.detail),
                 "reads_shown": _n(sum(c.reads or 0 for c in sel.shown), 1),
                 "summary": rt.summary_sentence(counts, sel.floor)},
    }


def _footer():
    return html.escape(rt.STRINGS["ui.footer"].format(
        present=f"{rt.DISPLAY_PRESENT:g}", weak=f"{rt.DISPLAY_WEAK:g}",
        p10=f"{rt.DISPLAY_10A_PRESENT:g}", w10=f"{rt.DISPLAY_10A_WEAK:g}",
        psize=f"{rt.DISPLAY_SIZE_PRESENT:g}", wsize=f"{rt.DISPLAY_SIZE_WEAK:g}"))


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
