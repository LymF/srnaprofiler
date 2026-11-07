"""Model selection (spec 7.3) from folds.jsonl; writes a markdown comparison table."""
from collections import defaultdict
from pathlib import Path

import numpy as np

from benchmark.ml.cv import read_lines
from benchmark.ml.stats import corrected_resampled_ttest

SIMPLICITY = ["logreg", "kmeans", "gmm", "hier", "hdbscan", "labelprop", "labelspread",
              "svm", "rf", "selftrain", "lgbm", "xgb"]
TOL = 0.01
ALPHA = 0.05
SCALARS = ["macro_f1", "ambiguous_share", "real_repro_gap", "control_real_rate",
           "rules_control_rate", "ece"]


class SelectionError(Exception):
    pass


def _get(rec, group, name, cls=None):
    m = rec.get(group) or {}
    if cls is None:
        return m.get(name)
    return ((m.get("per_class") or {}).get(cls) or {}).get(name)


def _mean(xs):
    xs = [x for x in xs if x is not None]
    return float(np.mean(xs)) if xs else None


def _fold_series(recs, units, fn):
    """Per-unit seed-averaged values of fn(rec) over `units` (None where unavailable)."""
    by = defaultdict(list)
    for r in recs:
        if r["unit"] in units:
            by[r["unit"]].append(fn(r))
    return {u: _mean(v) for u, v in by.items()}


def _paired(recs_a, recs_b, fn):
    """Aligned per-fold seed-averaged values of two record sets on their shared units."""
    units = {r["unit"] for r in recs_a} & {r["unit"] for r in recs_b}
    sa, sb = _fold_series(recs_a, units, fn), _fold_series(recs_b, units, fn)
    keep = sorted(u for u in units if sa.get(u) is not None and sb.get(u) is not None)
    return [sa[u] for u in keep], [sb[u] for u in keep], keep


def _sizes(recs, units):
    rs = [r for r in recs if r["unit"] in units]
    return (float(np.mean([r["n_train"] for r in rs])), float(np.mean([r["n_val"] for r in rs])))


def pvalue(recs_a, recs_b, metric="macro_f1"):
    recs_a = [r for r in recs_a if not is_xspecies(r["unit"])]
    recs_b = [r for r in recs_b if not is_xspecies(r["unit"])]
    a, b, units = _paired(recs_a, recs_b, lambda r: _get(r, "metrics", metric))
    if len(a) < 2:
        return float("nan")
    n_tr, n_te = _sizes(recs_a, set(units))
    return corrected_resampled_ttest(a, b, n_tr, n_te)[1]


def is_xspecies(unit):
    return str(unit).startswith("xspecies:")


def _aligned(recs, rules, fn, fn_rules=None):
    """(alg_values, rules_values) per unit where both seed-averaged values exist."""
    units = {r["unit"] for r in recs} & {r["unit"] for r in rules}
    sa = _fold_series(recs, units, fn)
    sb = _fold_series(rules, units, fn_rules or fn)
    keep = sorted(u for u in units if sa.get(u) is not None and sb.get(u) is not None)
    return [sa[u] for u in keep], [sb[u] for u in keep]


def eligibility(recs, rules):
    """(eligible, [failed criteria]) of one algorithm's main-variant records vs the rules'.
    Cross-species units are ignored. Every comparison uses exactly the folds where both sides
    have a value (seed-averaged per fold)."""
    recs = [r for r in recs if not is_xspecies(r["unit"])]
    rules = [r for r in rules if not is_xspecies(r["unit"])]
    why = []
    classes = sorted({c for r in rules for c in (r["metrics"].get("per_class") or {})})
    for c in classes:
        for name in ("precision", "recall"):
            fn = lambda r, c=c, name=name: _get(r, "metrics", name, c)
            mine, ref = _aligned(recs, rules, fn)
            if not ref:
                if name == "recall" and any(fn(r) is not None for r in rules):
                    why.append(f"{c} {name} missing")
                continue          # class counts only where rules' value exists
            if np.mean(mine) < np.mean(ref) - TOL:
                why.append(f"{c} {name} {np.mean(mine):.3f} < rules {np.mean(ref):.3f} - {TOL}")
    g = lambda n: (lambda r: _get(r, "metrics", n))
    amb, ramb = _aligned(recs, rules, g("ambiguous_share"))
    if not amb or not np.mean(amb) < np.mean(ramb):
        why.append(f"ambiguous_share {_mean(amb)} not < rules {_mean(ramb)}")
    gap = _mean([v for v in _fold_series(recs, {r["unit"] for r in recs},
                                         g("real_repro_gap")).values()])
    if gap is None or not gap > 0:
        why.append(f"real_repro_gap {gap} not > 0")
    ctl, rctl = _aligned(recs, recs, g("control_real_rate"), g("rules_control_rate"))
    if not ctl or not np.mean(ctl) < np.mean(rctl):
        why.append(f"control_real_rate {_mean(ctl)} not < rules_control_rate {_mean(rctl)}")
    return not why, why


def pick_winner(eligible: dict, rules):
    """Highest mean macro_f1 on the folds shared with rules; ties within p >= 0.05 against the
    runner-up go to the simpler algorithm."""
    if not eligible:
        return None, None, None
    f1 = {}
    for a, rs in eligible.items():
        mine, _ = _aligned([r for r in rs if not is_xspecies(r["unit"])],
                           [r for r in rules if not is_xspecies(r["unit"])],
                           lambda r: _get(r, "metrics", "macro_f1"))
        f1[a] = _mean(mine) or 0.0
    ranked = sorted(f1, key=lambda a: -f1[a])
    win = ranked[0]
    if len(ranked) == 1:
        return win, None, None
    run = ranked[1]
    p = pvalue(eligible[win], eligible[run])
    if (np.isnan(p) or p >= ALPHA) and SIMPLICITY.index(run) < SIMPLICITY.index(win):
        win = run
    return win, run, p


def _ms(xs, nd=3):
    xs = [x for x in xs if x is not None]
    if not xs:
        return "n/a"
    return f"{np.mean(xs):.{nd}f} ± {np.std(xs, ddof=1) if len(xs) > 1 else 0.0:.{nd}f}"


def _table(groups, rules, elig, group, header_extra, classes, xs=False):
    cols = SCALARS
    head = ["algorithm", "variant", "n"] + cols
    if group == "metrics":
        head += [f"{c} {m}" for c in classes for m in ("precision", "recall")]
        head += ["p vs rules", "eligible"]
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for (variant, alg), recs in sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        recs = [r for r in recs if is_xspecies(r["unit"]) == xs]
        if not recs:
            continue
        row = [alg, variant, str(len(recs))]
        row += [_ms([_get(r, group, s) for r in recs]) for s in cols]
        if group == "metrics":
            row += [_ms([_get(r, "metrics", m, c) for r in recs])
                    for c in classes for m in ("precision", "recall")]
            if xs:
                row += ["-", "excluded"]
            elif alg == "rules":
                row += ["-", "baseline"]
            else:
                p = pvalue(recs, rules)
                e = elig.get((variant, alg))
                row += [f"{p:.3g}", "-" if e is None else ("yes" if e[0] else "no")]
        lines.append("| " + " | ".join(row) + " |")
    return lines


def failed_combos(lines, algorithms=None):
    """Combos with an error line and no later success line."""
    ok = {(r["profile"], r["variant"], r["algorithm"], int(r["seed"]), r["unit"])
          for r in lines if "error" not in r}
    bad = {}
    for r in lines:
        k = (r["profile"], r["variant"], r["algorithm"], int(r["seed"]), r["unit"])
        if "error" in r and k not in ok and (algorithms is None or r["algorithm"] in algorithms
                                             or r["algorithm"] == "rules"):
            bad[k] = r["error"]
    return bad


def missing_combos(lines, algorithms=None, seeds=None):
    """Lines (strings) for every ranked (profile, variant, algorithm) lacking a success line
    for some (seed, unit): seeds must be exactly {0..max} with max+1 >= the seeds of the run,
    units those of the profile's rules baseline."""
    ok = [r for r in lines if "error" not in r]
    need = int(seeds) if seeds else max([int(r.get("seeds", 1)) for r in ok] or [1])
    out = []
    for prof in sorted({r["profile"] for r in ok}):
        pr = [r for r in ok if r["profile"] == prof]
        units = sorted({r["unit"] for r in pr if r["algorithm"] == "rules"})
        for v, a in sorted({(r["variant"], r["algorithm"]) for r in pr
                            if a_ranked(r["algorithm"], algorithms)}):
            have = {(int(r["seed"]), r["unit"]) for r in pr
                    if r["variant"] == v and r["algorithm"] == a}
            seeds_present = {s for s, _ in have}
            want = range(1 if a == "rules" else max(need, max(seeds_present) + 1))
            for sd in want:
                for u in units:
                    if (sd, u) not in have:
                        out.append(f"{prof}/{v}/{a}: seed {sd} unit {u}")
    return out


def a_ranked(alg, algorithms):
    return alg == "rules" or algorithms is None or alg in algorithms


def select(folds: Path, out: Path, algorithms: list[str] | None = None,
           seeds: int | None = None) -> str:
    lines = read_lines(folds)
    bad = failed_combos(lines, algorithms)
    if bad:
        msg = "\n".join(f"  {k}: {e}" for k, e in sorted(bad.items()))
        raise SelectionError(f"{len(bad)} failed fold(s); rerun ml-cv to retry:\n{msg}")
    miss = missing_combos(lines, algorithms, seeds)
    if miss:
        raise SelectionError(f"{len(miss)} missing fold result(s); rerun ml-cv:\n  "
                             + "\n  ".join(miss[:50]))
    lines = [r for r in lines if "error" not in r
             and (algorithms is None or r["algorithm"] in algorithms or r["algorithm"] == "rules")]
    out_md = ["# ML classifier selection", ""]
    for profile in sorted({r["profile"] for r in lines}):
        pl = [r for r in lines if r["profile"] == profile]
        groups = defaultdict(list)
        for r in pl:
            groups[(r["variant"], r["algorithm"])].append(r)
        rules = groups.get(("main", "rules"))
        if not rules:
            raise SelectionError(f"profile {profile}: no rules baseline lines in {folds}")
        elig = {(v, a): eligibility(rs, rules) for (v, a), rs in groups.items()
                if v == "main" and a != "rules"}
        pool = {a: groups[("main", a)] for (v, a), e in elig.items() if e[0]}
        win, run, p = pick_winner(pool, rules)
        classes = sorted({c for r in rules for c in (r["metrics"].get("per_class") or {})})
        out_md += [f"## {profile}", "",
                   "Values: mean ± sd over folds×seeds. p: corrected resampled t-test of "
                   "macro_f1 against rules (per-fold, seed-averaged).", ""]
        out_md += _table(groups, rules, elig, "metrics", None, classes)
        out_md += ["", "### complete_rules (rule class kept where not ambiguous, ML elsewhere; "
                   "not used for eligibility)", ""]
        out_md += _table(groups, rules, elig, "complete_rules_metrics", None, classes)
        if any(is_xspecies(r["unit"]) for r in pl):
            out_md += ["", "### Cross-species folds (own group; excluded from eligibility, "
                       "ranking and t-tests)", ""]
            out_md += _table(groups, rules, elig, "metrics", None, classes, xs=True)
        out_md += ["", "### Eligibility (variant main)", ""]
        for (v, a), (ok, why) in sorted(elig.items(), key=lambda kv: kv[0][1]):
            out_md.append(f"- {a}: " + ("eligible" if ok else "not eligible: " + "; ".join(why)))
        out_md.append("")
        if win is None:
            out_md.append("**Winner: none (no eligible algorithm)**")
        else:
            extra = "" if run is None else f" (runner-up {run}, p = {p:.3g})"
            out_md.append(f"**Winner: {win}**{extra}")
        out_md.append("")
    text = "\n".join(out_md)
    Path(out).parent.mkdir(parents=True, exist_ok=True)
    Path(out).write_text(text)
    return text
