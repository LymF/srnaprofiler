"""Model selection (spec 7.3) from folds.jsonl; writes a markdown comparison table."""
from collections import defaultdict
from pathlib import Path

import numpy as np

from benchmark.ml.cv import read_lines
from benchmark.ml.stats import corrected_resampled_ttest

SIMPLICITY = ["logreg", "kmeans", "gmm", "hier", "hdbscan", "labelprop", "labelspread",
              "svm", "rf", "selftrain", "lgbm", "xgb"]
CRITERION = "v2"
SCALARS = ["macro_f1", "ambiguous_share", "real_repro_gap", "control_real_rate",
           "rules_control_rate", "unann_supported_rate", "rules_unann_supported_rate", "ece"]


def selection_config() -> dict:
    """{class_floor, alpha} of criterion v2 from ml/config.yaml."""
    from benchmark.ml.dataset import load_ml_config
    sel = load_ml_config()["selection"]
    return {"class_floor": float(sel["class_floor"]), "alpha": float(sel["alpha"])}


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


def eligibility(recs, rules, class_floor=None, alpha=None):
    """(eligible, [failed criteria]) of one algorithm's main-variant records vs the rules'
    (criterion v2, spec section 3):
      1. macro_f1 mean > rules' and corrected resampled t-test p < alpha;
      2. no class loses more than class_floor precision or recall;
      3. unann_supported_rate mean >= the rules' (independent biology of unannotated calls);
      4. control_real_rate mean < rules_control_rate mean;
      5. ambiguous_share mean < rules'.
    real_repro_gap is reported only. Cross-species units are ignored. Every comparison uses
    exactly the folds where both sides have a value (seed-averaged per fold)."""
    cfg = selection_config()
    class_floor = cfg["class_floor"] if class_floor is None else class_floor
    alpha = cfg["alpha"] if alpha is None else alpha
    recs = [r for r in recs if not is_xspecies(r["unit"])]
    rules = [r for r in rules if not is_xspecies(r["unit"])]
    why = []
    g = lambda n: (lambda r: _get(r, "metrics", n))
    f1, rf1 = _aligned(recs, rules, g("macro_f1"))
    p = pvalue(recs, rules)
    if not f1 or not np.mean(f1) > np.mean(rf1) or np.isnan(p) or not p < alpha:
        why.append(f"macro_f1 {_mean(f1)} vs rules {_mean(rf1)}, p {p:.3g} (need higher and p < {alpha})")
    classes = sorted({c for r in rules for c in (r["metrics"].get("per_class") or {})})
    for c in classes:
        for name in ("precision", "recall"):
            fn = lambda r, c=c, name=name: _get(r, "metrics", name, c)
            mine, ref = _aligned(recs, rules, fn)
            if not ref:
                if name == "recall" and any(fn(r) is not None for r in rules):
                    why.append(f"{c} {name} missing")
                continue          # class counts only where rules' value exists
            if np.mean(mine) < np.mean(ref) - class_floor:
                why.append(f"{c} {name} {np.mean(mine):.3f} < rules {np.mean(ref):.3f} - {class_floor}")
    sup, rsup = _aligned(recs, recs, g("unann_supported_rate"), g("rules_unann_supported_rate"))
    if not sup or not np.mean(sup) >= np.mean(rsup):
        why.append(f"unann_supported_rate {_mean(sup)} not >= rules {_mean(rsup)}")
    ctl, rctl = _aligned(recs, recs, g("control_real_rate"), g("rules_control_rate"))
    if not ctl or not np.mean(ctl) < np.mean(rctl):
        why.append(f"control_real_rate {_mean(ctl)} not < rules_control_rate {_mean(rctl)}")
    amb, ramb = _aligned(recs, rules, g("ambiguous_share"))
    if not amb or not np.mean(amb) < np.mean(ramb):
        why.append(f"ambiguous_share {_mean(amb)} not < rules {_mean(ramb)}")
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
    if (np.isnan(p) or p >= selection_config()["alpha"]) and SIMPLICITY.index(run) < SIMPLICITY.index(win):
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


def missing_combos(lines, algorithms=None, seeds=None, variants=None):
    """Messages for every requested (profile, variant, algorithm) lacking a success line for
    some (seed, unit): seeds must be {0..max} with max+1 >= the run's seeds; units are all
    units seen in the profile (so rules lacking a unit any algorithm has is flagged too).
    Requested algorithms default to every non-rules entry of REGISTRY, variants to
    main,stack,repro (rules only exist as main, seed 0)."""
    from benchmark.ml.cv import REGISTRY, VARIANTS
    ok = [r for r in lines if "error" not in r]
    need = int(seeds) if seeds else max([int(r.get("seeds", 1)) for r in ok] or [1])
    algs = [a for a in (algorithms or sorted(REGISTRY)) if a != "rules"]
    variants = list(variants or VARIANTS)
    out = []
    for prof in sorted({r["profile"] for r in ok}):
        pr = [r for r in ok if r["profile"] == prof]
        units = sorted({r["unit"] for r in pr})
        combos = [("main", "rules")] + [(v, a) for v in variants for a in algs]
        for v, a in combos:
            have = {(int(r["seed"]), r["unit"]) for r in pr
                    if r["variant"] == v and r["algorithm"] == a}
            if not have:
                out.append(f"{prof}/{v}/{a}: no results")
                continue
            top = 1 if a == "rules" else max(need, max(s for s, _ in have) + 1)
            out += [f"{prof}/{v}/{a}: seed {sd} unit {u}"
                    for sd in range(top) for u in units if (sd, u) not in have]
    return out


def select(folds: Path, out: Path, algorithms: list[str] | None = None,
           seeds: int | None = None, variants: list[str] | None = None) -> str:
    lines = read_lines(folds)
    bad = failed_combos(lines, algorithms)
    if bad:
        msg = "\n".join(f"  {k}: {e}" for k, e in sorted(bad.items()))
        raise SelectionError(f"{len(bad)} failed fold(s); rerun ml-cv to retry:\n{msg}")
    miss = missing_combos(lines, algorithms, seeds, variants)
    if miss:
        raise SelectionError(f"{len(miss)} missing fold result(s); rerun ml-cv:\n  "
                             + "\n  ".join(miss[:50]))
    lines = [r for r in lines if "error" not in r
             and (algorithms is None or r["algorithm"] in algorithms or r["algorithm"] == "rules")]
    cfg = selection_config()
    out_md = ["# ML classifier selection", "",
              f"Selection criterion {CRITERION}: class_floor = {cfg['class_floor']}, "
              f"alpha = {cfg['alpha']}. Eligible = macro-F1 above rules with p < alpha; no class "
              "loses more than class_floor precision or recall; supported share of unannotated "
              "calls (repro_frac >= 0.5, plus TE overlap for piRNA where a TE file is "
              "configured) >= rules'; control rate below rules'; ambiguous share below rules'. "
              "real_repro_gap is reported, not gated.", ""]
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
