"""Class checklist for the cluster page ("Why this class").

One row per check the rules apply to the class (rules.yaml / rules.plant.yaml, as read by
`report_data.class_specs`): the size window, each defining evidence and the optional cluster
criteria. The decider weighs the evidence in a score and has no cutoff on it; what it does
gate hard is the size window, the cluster criteria and `phase_max_pvalue`. An evidence row
therefore passes when its effect strength (`e_*`, a 95% Wilson lower bound against the
signature's null, see evidence.effect_strength) is above 0, i.e. above chance.

Payload per cluster: [class, [status, observed text], ...], rows in `skeleton(kingdom)[class]`
order. `skeleton` ([[name key, rule text], ...] per class) goes into meta once.
"""
import functools

import report_data
import report_text as rt

PASS, FAIL, NA, UNTESTED = 1, 0, 2, 3
_NAME = {"1u": "1u", "5a": "5a", "fractop_bid": "strand", "pingpong": "pingpong",
         "phase": "phasing", "size": "size"}
_E = ("1u", "5a", "fractop_bid", "pingpong", "phase")


def _g(x):
    return f"{x:g}"


@functools.lru_cache(maxsize=None)
def _plan(kingdom):
    """{class: [row dict]}: size first, then the other slots, then the cluster criteria."""
    out = {}
    for cls, spec in report_data.class_specs(kingdom).items():
        rows = [{"kind": "size", "names": ["size"]}] if any("size" in s for s in spec["slots"]) else []
        rows += [{"kind": "ev", "names": s} for s in spec["slots"] if "size" not in s]
        if spec["min_length"] is not None:
            rows.append({"kind": "length"})
        if spec["density_p"] is not None:
            rows.append({"kind": "density"})
        out[cls] = rows
    return out


def _key(row, spec):
    if row["kind"] != "ev":
        return row["kind"]
    key = "_".join(_NAME[n] for n in row["names"])
    return key + "_p" if key == "phasing" and spec["phase_max_p"] is not None else key


def _rule(row, spec):
    S = rt.STRINGS
    k = row["kind"]
    if k == "size":
        lo, hi = spec["size"]
        return S["ck.thr.size1" if lo == hi else "ck.thr.size"].format(lo=lo, hi=hi)
    if k == "length":
        return S["ck.thr.length"].format(n=_g(spec["min_length"]))
    if k == "density":
        return S["ck.thr.density"].format(p=_g(spec["density_p"]))
    return S["ck.thr." + _key(row, spec)].format(p=_g(spec["phase_max_p"] or 0))


def _name_key(row):
    k = row["kind"]
    key = "size" if k == "size" else k if k != "ev" else "_".join(_NAME[n] for n in row["names"])
    return f"check.{key}.title"


@functools.lru_cache(maxsize=None)
def skeleton(kingdom):
    """{class: [[name key in STRINGS, rule text], ...]} for every class with a checklist."""
    specs = report_data.class_specs(kingdom)
    return {cls: [[_name_key(r), _rule(r, specs[cls])] for r in rows]
            for cls, rows in _plan(kingdom).items()}


def _strength(c, spec, name):
    """e_* of one evidence as the decider reads it; None when the table lacks it."""
    if name == "phase" and spec["phase_max_p"] is not None:
        p = c.phase_pvalue
        if p is None or p != p or p > spec["phase_max_p"]:
            return 0.0
    return (c.e or {}).get(name)


def _obs(c, name):
    """Observed value of one evidence, worded (strength appended when the table has it)."""
    S = rt.STRINGS
    na = S["ui.na"]
    e = (c.e or {}).get(name)
    if name == "1u":
        txt = None if c.pct_1u is None else S["ck.obs.1u"].format(pct=round(c.pct_1u * 100))
    elif name == "5a":
        txt = None if c.pct_5a is None else S["ck.obs.5a"].format(pct=round(c.pct_5a * 100))
    elif name == "fractop_bid":
        txt = None if c.frac_top is None else S["ck.obs.strand"].format(
            pct=round(min(c.frac_top, 1 - c.frac_top) * 100))
    elif name == "pingpong":
        txt = None if c.pingpong_z is None else S["ck.obs.pingpong"].format(z=f"{c.pingpong_z:.1f}")
    else:  # phase
        if c.phase_pvalue is None:
            txt = None
        else:
            n = "" if c.phase_number is None else round(c.phase_number)
            txt = S["ck.obs.phasing"].format(p=f"{c.phase_pvalue:.2g}", n=n)
    if txt is None:
        return na
    return txt + (S["ck.obs.strength"].format(e=f"{e:.2f}") if e is not None else "")


def _ev_status(c, spec, names):
    vals = [_strength(c, spec, n) for n in names]
    if any(v is not None and v > 0 for v in vals):
        return PASS
    return NA if any(v is None for v in vals) else FAIL


def checklist(c, kingdom="animal"):
    """[class, [status, observed], ...] for a cluster, or None when there is nothing to list.

    A resolved class lists its own checks; an ambiguous cluster the checks of the nearest
    class named in its evidence trail (none for an ML call without one). miRNA comes from
    the miRDeep2 hairpin search and has no read-signature checks.
    """
    S = rt.STRINGS
    amb = c.cls == "ambiguous"
    cls = rt.trail_nearest(c) if amb else c.cls
    specs = report_data.class_specs(kingdom)
    if cls not in specs:
        return None
    spec = specs[cls]
    failed = set()
    m = rt._CRITERION.search(c.evidence or "") if amb else None
    if m:
        failed = {n.strip() for n in m.group(1).split(",")}
    z = c.dicer_call
    in_window = z is not None and spec["size"] is not None and spec["size"][0] <= z <= spec["size"][1]
    out = [cls]
    for row in _plan(kingdom)[cls]:
        k = row["kind"]
        if k == "size":
            out.append([PASS if in_window else FAIL,
                        S["ck.obs.nosize"] if z is None else S["ck.obs.size"].format(z=z)])
        elif k == "ev":
            obs = " / ".join(_obs(c, n) for n in row["names"])
            out.append([_ev_status(c, spec, row["names"]), obs])
        else:
            # criteria are only tested for a locus whose size fits the window
            if amb and not in_window:
                st = UNTESTED
            else:
                st = FAIL if k in failed else PASS
            obs = S["ck.obs.span"].format(n=_g(c.end - c.start)) \
                if k == "length" and c.start is not None and c.end is not None else S["ui.na"]
            out.append([st, obs])
    return out
