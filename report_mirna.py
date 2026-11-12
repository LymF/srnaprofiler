"""miRNA hairpin data for the cluster page of the HTML report.

Reads `<prefix>.mirna.tsv` (see mirna_report.MIRNA_COLS), matches its rows to the clusters
classified miRNA, and builds the per-page data: the precursor with its mature, loop and star
segments, the fold (RNAfold when installed), the reads inside the precursor, the numbers and
the plant criteria of Axtell and Meyers 2018.

Coordinates: start/end are 0-based half-open like clusters.tsv. The precursor sequence is the
one miRScore was given: the engine's hairpin re-cut from the genome with FLANK nt on each side
(clipped at a chromosome end), written in the strand's orientation, so for a minus-strand
precursor position 0 is the highest genomic coordinate. mature/star start-end are 1-based
inclusive positions in that sequence; they are checked against the sequences and corrected
when they do not fit.
"""
import csv
import importlib
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass

import report_text as rt

FLANK = 25          # mirna.hairpin_with_flank's default
READ_ROWS = 25      # read variants drawn; the rest are counted
# Axtell and Meyers 2018 (Plant Cell 30:272-284), the three that the table can answer.
PLANT_MAX_PRECURSOR = 300
PLANT_MAX_MISMATCHES = 5
PLANT_MIN_PRECISION = 75.0

_REQUIRED = ("name", "chrom", "start", "end", "strand")
_NA = ("", "NA", "N", "nan")


def _s(x):
    return None if x is None or x.strip() in _NA else x.strip()


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def _i(x):
    v = _f(x)
    return None if v is None else int(v)


def _rna(seq):
    return None if seq is None else seq.upper().replace("T", "U")


@dataclass
class Row:
    name: str
    chrom: str
    start: int
    end: int
    strand: str
    mature_seq: str = None
    star_seq: str = None
    precursor_seq: str = None
    precursor_len: int = None
    mature_start: int = None
    mature_end: int = None
    star_start: int = None
    star_end: int = None
    result: str = None
    flags: str = None
    mismatches: int = None
    bulges: int = None
    mature_reads: float = None
    star_reads: float = None
    total_reads: float = None
    precision: float = None
    mirdeep_score: str = None
    mirdeep_prob: str = None
    rfam_alert: str = None
    randfold: str = None
    mirbase_hit: str = None
    mirbase_same_seed: str = None


def load(path):
    """Rows of a mirna.tsv. A file without the locus columns is an error (ValueError)."""
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in _REQUIRED if c not in (rd.fieldnames or [])]
        if missing:
            raise ValueError(f"{path} lacks the column(s) {', '.join(missing)}; "
                             f"is it a <prefix>.mirna.tsv?")
        out = []
        for r in rd:
            if _i(r.get("start")) is None or _i(r.get("end")) is None:
                continue
            out.append(Row(
                name=r["name"], chrom=r["chrom"], start=_i(r["start"]), end=_i(r["end"]),
                strand=(r.get("strand") or ".").strip(),
                mature_seq=_rna(_s(r.get("mature_seq"))), star_seq=_rna(_s(r.get("star_seq"))),
                precursor_seq=_rna(_s(r.get("precursor_seq"))),
                precursor_len=_i(r.get("precursor_len")),
                mature_start=_i(r.get("mature_start")), mature_end=_i(r.get("mature_end")),
                star_start=_i(r.get("star_start")), star_end=_i(r.get("star_end")),
                result=_s(r.get("result")), flags=_s(r.get("flags")),
                mismatches=_i(r.get("mismatches")), bulges=_i(r.get("bulges")),
                mature_reads=_f(r.get("mature_reads")), star_reads=_f(r.get("star_reads")),
                total_reads=_f(r.get("total_reads")), precision=_f(r.get("precision")),
                mirdeep_score=_s(r.get("mirdeep_score")), mirdeep_prob=_s(r.get("mirdeep_prob")),
                rfam_alert=_s(r.get("rfam_alert")), randfold=_s(r.get("randfold")),
                mirbase_hit=_s(r.get("mirbase_hit")),
                mirbase_same_seed=_s(r.get("mirbase_same_seed"))))
        return out


def match(clusters, rows, min_fraction=0.5):
    """{cluster name: Row} for the clusters classified miRNA.

    The rule of `classify --mirna-loci` (classify.mark_mirna): same chromosome and the
    miRNA interval covers at least `min_fraction` of the cluster's span. With several
    candidates the largest overlap wins, a Pass before a Fail, then the table order.
    """
    by_chrom = {}
    for i, r in enumerate(rows):
        by_chrom.setdefault(r.chrom, []).append((i, r))
    out = {}
    for c in clusters:
        if c.cls != "miRNA" or c.start is None or c.end is None or c.end <= c.start:
            continue
        span = c.end - c.start
        best = None
        for i, r in by_chrom.get(c.chrom, []):
            ov = min(c.end, r.end) - max(c.start, r.start)
            if ov <= 0 or ov / span < min_fraction:
                continue
            key = (-ov, 0 if r.result == "Pass" else 1, i)
            if best is None or key < best[0]:
                best = (key, r)
        if best:
            out[c.name] = best[1]
    return out


# ---- geometry ----

def window(row):
    """(lo, hi) genomic half-open window of the precursor sequence, or None.

    The sequence is the hairpin plus up to FLANK nt each side; clipping happens only at a
    chromosome end, so when the two flanks are not both FLANK the shorter one is the left
    one if the hairpin starts within FLANK of position 0 (else the right one).
    """
    seq = row.precursor_seq
    if not seq:
        return None
    L = row.precursor_len if row.precursor_len else len(seq)
    if L != len(seq):
        return None
    extra = L - (row.end - row.start)
    if extra < 0 or extra > 2 * FLANK:
        return None
    left = FLANK if extra == 2 * FLANK else min(FLANK, row.start)
    right = extra - left
    if right < 0 or right > FLANK:
        return None
    return row.start - left, row.end + right


def to_precursor(ref_start, ref_end, strand, lo, hi):
    """[a, b) in precursor coordinates of the genomic interval [ref_start, ref_end)."""
    if strand == "-":
        return hi - ref_end, hi - ref_start
    return ref_start - lo, ref_end - lo


def _locate(seq, sub, start, end):
    """0-based half-open (a, b) of `sub` in `seq`: the 1-based inclusive start/end given,
    else the same read as 0-based, else the first occurrence; None when not found."""
    if not sub:
        return None
    n = len(sub)
    if start is not None and end is not None:
        for a in (start - 1, start):
            if a >= 0 and seq[a:a + n] == sub and a + n in (end, end + 1):
                return a, a + n
    i = seq.find(sub)
    return None if i < 0 else (i, i + n)


def segments(row):
    """{"m": (a, b) | None, "s": ..., "loop": ...}: 0-based half-open in the precursor.

    The loop is the stretch between the two segments (needs both).
    """
    seq = row.precursor_seq
    if not seq:
        return {"m": None, "s": None, "loop": None}
    m = _locate(seq, row.mature_seq, row.mature_start, row.mature_end)
    s = _locate(seq, row.star_seq, row.star_start, row.star_end)
    loop = None
    if m and s:
        lo, hi = (m[1], s[0]) if m[0] <= s[0] else (s[1], m[0])
        loop = (lo, hi) if hi > lo else None
    return {"m": m, "s": s, "loop": loop}


# ---- structure ----

def _import_rna():
    try:
        return importlib.import_module("RNA")
    except ImportError:
        return None


def fold(seq):
    """(dot-bracket, mfe) from RNAfold (the executable, else the python module).

    Returns (None, None) when neither is available or the fold fails. Nothing is cached or
    left on disk.
    """
    if not seq:
        return None, None
    exe = shutil.which("RNAfold")
    if exe:
        try:
            with tempfile.TemporaryDirectory() as tmp:
                p = subprocess.run([exe, "--noPS"], input=seq + "\n", capture_output=True,
                                   text=True, timeout=120, cwd=tmp)
            lines = p.stdout.strip().splitlines()
            if p.returncode == 0 and len(lines) >= 2:
                m = re.match(r"^([().]+)\s*\(\s*(-?\d+(?:\.\d+)?)\s*\)\s*$", lines[1].strip())
                if m and len(m.group(1)) == len(seq):
                    return m.group(1), float(m.group(2))
        except (OSError, subprocess.SubprocessError, ValueError):
            pass
        return None, None
    mod = _import_rna()
    if mod is not None:
        try:
            db, mfe = mod.fold(seq)
            if len(db) == len(seq):
                return db, float(mfe)
        except Exception:
            pass
    return None, None


def fold_tool_missing():
    """True when neither RNAfold nor the python module is there."""
    return shutil.which("RNAfold") is None and _import_rna() is None


# ---- numbers and criteria ----

def numbers(row):
    """[[key, value]] for the numbers table, in display order; None when not in the table."""
    core = row.end - row.start
    return [["result", row.result], ["mature_reads", row.mature_reads],
            ["star_reads", row.star_reads], ["precision", row.precision],
            ["mismatches", row.mismatches], ["bulges", row.bulges],
            ["precursor_len", row.precursor_len], ["locus_len", core],
            ["mirdeep_score", row.mirdeep_score], ["mirdeep_prob", row.mirdeep_prob],
            ["randfold", row.randfold], ["rfam_alert", row.rfam_alert],
            ["mirbase_hit", row.mirbase_hit], ["mirbase_same_seed", row.mirbase_same_seed]]


def criteria(row):
    """[[key, status, observed text, rule text]], status 1 met, 0 not met, 2 not available.

    Only the criteria the table can answer. The precursor length is the engine's hairpin
    (end - start), not the sequence drawn, which carries the flanks.
    """
    S = rt.STRINGS
    na = S["ui.na"]
    core = row.end - row.start
    out = [["len", 1 if core <= PLANT_MAX_PRECURSOR else 0,
            S["mi.obs.nt"].format(n=core), S["mi.rule.len"].format(n=PLANT_MAX_PRECURSOR)]]
    mm = row.mismatches
    out.append(["mm", 2 if mm is None else 1 if mm <= PLANT_MAX_MISMATCHES else 0,
                na if mm is None else str(mm), S["mi.rule.mm"].format(n=PLANT_MAX_MISMATCHES)])
    p = row.precision
    out.append(["prec", 2 if p is None else 1 if p >= PLANT_MIN_PRECISION else 0,
                na if p is None else S["mi.obs.pct"].format(p=f"{p:g}"),
                S["mi.rule.prec"].format(p=f"{PLANT_MIN_PRECISION:g}")])
    return out


# ---- payload ----

def hairpin_spec(row):
    """What the BAM pass needs: (lo, hi, strand, precursor, mature (a, b)|None, star|None),
    plain tuples so it crosses a process boundary; None when the precursor cannot be placed."""
    win = window(row)
    if win is None:
        return None
    seg = segments(row)
    return (win[0], win[1], "-" if row.strand == "-" else "+", row.precursor_seq,
            seg["m"], seg["s"])


def figure_data(row, folded, hp):
    """The `mi` entry of a detail page.

    folded: (dot-bracket|None, mfe|None, tool_missing bool); hp: the reads dict of
    report_signals.HairpinReads.build() or None.
    Keys: n name, st strand, L precursor length, seq (RNA) or absent, m/s [a, b] segments,
    lp [a, b] loop, db dot-bracket, mfe, ns why there is no structure (m missing, f failed),
    rd [[a, b, sequence (case = match), reads, kind m/s/""]], cov per-base depth, more/mr
    variants and reads not drawn, anti/ovh reads not shown, lo genomic window start,
    t numbers [[key, value]], pc criteria [[key, status, observed, rule]].
    """
    seg = segments(row)
    d = {"n": row.name, "st": "-" if row.strand == "-" else "+",
         "t": numbers(row), "pc": criteria(row)}
    if not row.precursor_seq:
        return d
    db, mfe, missing = folded if folded else (None, None, False)
    d.update({"L": len(row.precursor_seq), "seq": row.precursor_seq})
    for k, key in (("m", "m"), ("s", "s"), ("lp", "loop")):
        if seg[key]:
            d[k] = list(seg[key])
    if db:
        d["db"] = db
        if mfe is not None:
            d["mfe"] = mfe
    else:
        d["ns"] = "m" if missing else "f"
    win = window(row)
    if win is None:
        d["nowin"] = 1
    elif hp is not None:
        d.update({k: hp[k] for k in ("rd", "cov", "more", "mr", "anti", "ovh") if k in hp})
        d["lo"] = win[0]
    return d
