"""Library quality strip for the Overview tab: read-length flags and trimming numbers.

Flags in the miRTrace style, each a number and an ok / check status. The 20-25 nt share
and its 25% floor follow the brief; the short-read cutoff (more than 25% of reads under
18 nt) and the "dominant length inside the sizes the classifier keeps" test are display-only
choices of this report, not miRTrace values. Unmapped reads are not counted, so every
share is of mapped reads.
"""
import json
import os
import re

import report_signals
from utils import log

SHARE_LO, SHARE_HI, SHARE_MIN = 20, 25, 0.25
SHORT_BELOW, SHORT_MAX = 18, 0.25
KEPT_RANGE = (18, 35)  # extractor dicermin / dicermax of both shipped rules files
_LO, _HI = report_signals.LEN_LO, report_signals.LEN_HI
_NUM = r"([\d,]+)"
_BEFORE = re.compile(r"Reads before:\s*" + _NUM)
_AFTER = re.compile(r"Reads after:\s*" + _NUM)


def _label(slot):
    return "<%d" % _LO if slot == 0 else ">%d" % _HI if slot == _HI - _LO + 2 else slot + _LO - 1


def build(lengths, trim=None):
    """Payload of the strip, or None without reads.

    {"c": [27 weighted counts, 0.1 precision], "u": unmapped records, "fl": [[key, value,
    status]] (status 1 ok, 0 check; keys share, dominant, short), "tr": {b, a, ad}}.
    """
    if not lengths:
        return None
    c = lengths["counts"]
    total = sum(c)
    if total <= 0:
        return None

    def at(n):
        return c[n - _LO + 1]
    share = sum(at(n) for n in range(SHARE_LO, SHARE_HI + 1)) / total
    short = (c[0] + sum(at(n) for n in range(_LO, SHORT_BELOW))) / total
    top = max(range(len(c)), key=lambda i: (c[i], -i))  # ties: the shorter length
    dom = _label(top)
    dom_ok = isinstance(dom, int) and KEPT_RANGE[0] <= dom <= KEPT_RANGE[1]
    out = {"c": [round(v, 1) for v in c], "u": lengths.get("unmapped"),
           "fl": [["share", round(share, 4), int(share >= SHARE_MIN)],
                  ["dominant", dom, int(dom_ok)],
                  ["short", round(short, 4), int(short <= SHORT_MAX)]]}
    if trim:
        out["tr"] = trim
    return out


def _int(s):
    return int(s.replace(",", ""))


def load_trim(fastp=None, qc_summary=None):
    """{"b": reads before, "a": after, "ad": adapter-trimmed} from fastp.json, with
    before/after filled from qc_summary.txt when fastp lacks them; None when nothing is found.
    An unreadable file is skipped with a warning."""
    out = {}
    if fastp:
        try:
            with open(fastp) as fh:
                j = json.load(fh)
            s = j.get("summary") or {}
            for k, sec in (("b", "before_filtering"), ("a", "after_filtering")):
                v = (s.get(sec) or {}).get("total_reads")
                if isinstance(v, (int, float)):
                    out[k] = int(v)
            v = (j.get("adapter_cutting") or {}).get("adapter_trimmed_reads")
            if isinstance(v, (int, float)):
                out["ad"] = int(v)
        except (OSError, ValueError, AttributeError) as e:
            log(f"report: could not read {fastp} ({e}); trimming numbers left out", "WARN")
    if qc_summary and not ("b" in out and "a" in out):
        try:
            with open(qc_summary) as fh:
                txt = fh.read()
            for k, rx in (("b", _BEFORE), ("a", _AFTER)):
                m = rx.search(txt)
                if m and k not in out:
                    out[k] = _int(m.group(1))
        except (OSError, ValueError) as e:
            log(f"report: could not read {qc_summary} ({e}); trimming numbers left out", "WARN")
    return {k: out[k] for k in ("b", "a", "ad") if k in out} or None


def trim_paths(clusters, fastp, qc_summary):
    """Explicit paths, else <prefix>.fastp.json / <prefix>.qc_summary.txt beside --clusters."""
    base = clusters[:-len(".clusters.tsv")] if clusters.endswith(".clusters.tsv") else None

    def pick(given, suffix):
        if given:
            return given
        cand = base + suffix if base else None
        return cand if cand and os.path.isfile(cand) else None
    return pick(fastp, ".fastp.json"), pick(qc_summary, ".qc_summary.txt")
