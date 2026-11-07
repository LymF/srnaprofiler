import json
import math

import numpy as np
import pandas as pd
import pytest
from scipy import stats

pytest.importorskip("sklearn")

from benchmark import __main__ as cli
from benchmark.ml import cv
from benchmark.ml.select import select, SelectionError
from benchmark.ml.stats import corrected_resampled_ttest
from mlfeatures import FEATURES

UNITS = [f"u{i}" for i in range(10)]


def test_corrected_ttest_known_value():
    z = np.arange(10, dtype=float)
    z = (z - z.mean()) / z.std(ddof=1)
    d = 0.01 + 0.005 * z                      # mean 0.01, sd 0.005, k = 10
    b = np.full(10, 0.5)
    t, p = corrected_resampled_ttest(list(b + d), list(b), n_train=90, n_test=10)
    expected_t = 0.01 / math.sqrt((1 / 10 + 10 / 90) * 0.005 ** 2)     # = 4.3529...
    assert expected_t == pytest.approx(4.35290, abs=1e-4)
    assert t == pytest.approx(expected_t)
    assert p == pytest.approx(2 * stats.t.sf(expected_t, 9))


def _m(f1, prec=0.9, rec=0.9, amb=0.1, gap=0.2, ctl=0.1, rctl=0.5):
    return {"per_class": {"siRNA": {"precision": prec, "recall": rec},
                          "piRNA": {"precision": 0.8, "recall": 0.8}},
            "macro_f1": f1, "ambiguous_share": amb, "real_repro_gap": gap,
            "control_real_rate": ctl, "rules_control_rate": rctl, "ece": 0.05}


def _line(alg, unit, seed, m, variant="main", **kw):
    return {"profile": "arthropod", "variant": variant, "algorithm": alg, "seed": seed,
            "unit": unit, "t": 0.5, "n_train": 90, "n_val": 10, "metrics": m,
            "complete_rules_metrics": m, **kw}


def _write(tmp_path, specs):
    """specs: alg -> (f1 list per unit, metric kwargs)."""
    rows = [_line("rules", u, 0, _m(0.5, amb=0.4, gap=None, ctl=None)) for u in UNITS]
    for alg, (f1s, kw) in specs.items():
        for u, f in zip(UNITS, f1s):
            rows.append(_line(alg, u, 0, _m(f, **kw)))
    path = tmp_path / "folds.jsonl"
    path.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    return path


def test_select_rejects_class_loss(tmp_path):
    f1 = [0.9 + 0.001 * i for i in range(10)]
    path = _write(tmp_path, {"xgb": (f1, {"prec": 0.88}),     # 0.02 below rules (0.9) on siRNA
                             "logreg": ([0.6 + 0.001 * i for i in range(10)], {})})
    text = select(path, tmp_path / "sel.md")
    assert "not eligible: siRNA precision" in text.split("xgb:")[1].split("\n")[0]
    assert "**Winner: logreg**" in text


def test_select_tie_prefers_simpler(tmp_path):
    rng = np.random.default_rng(1)
    noise = rng.normal(0, 0.02, 10)
    path = _write(tmp_path, {"rf": ([0.8 + x for x in noise], {}),
                             "xgb": ([0.8 + x for x in noise[::-1]], {})})
    text = select(path, tmp_path / "sel.md")
    assert "**Winner: rf**" in text


def test_select_clear_winner_beats_simplicity(tmp_path):
    f1 = [0.9 + 0.002 * i for i in range(10)]
    path = _write(tmp_path, {"xgb": (f1, {}), "rf": ([x - 0.2 + 0.003 * i for i, x in
                                                         enumerate(f1)], {})})
    assert "**Winner: xgb**" in select(path, tmp_path / "sel.md")


def test_select_refuses_on_error_lines(tmp_path, capsys):
    path = _write(tmp_path, {"rf": ([0.8] * 10, {})})
    with open(path, "a") as fh:
        fh.write(json.dumps({"profile": "arthropod", "variant": "main", "algorithm": "rf",
                             "seed": 1, "unit": "u3", "error": "ValueError('boom')"}) + "\n")
    out = tmp_path / "sel.md"
    assert cli.main(["ml-select", "--folds", str(path), "-o", str(out)]) == 1
    assert "boom" in capsys.readouterr().err and not out.exists()
    # an error for an algorithm not being ranked does not block
    assert cli.main(["ml-select", "--folds", str(path), "-o", str(out),
                     "--algorithms", "xgb"]) == 0


def _tiny_table(n=120, seed=0):
    rng = np.random.default_rng(seed)
    cls = rng.integers(0, 2, n)
    X = rng.normal(0, 0.5, (n, len(FEATURES))) + cls[:, None] * 3.0
    df = pd.DataFrame(X, columns=FEATURES)
    lab = np.array(["siRNA", "piRNA"], dtype=object)[cls]
    df["label"] = np.where(rng.random(n) < 0.5, lab, None)
    df["rule_cls"] = np.where(rng.random(n) < 0.3, "ambiguous", lab)
    df["carved"] = False
    df["is_control"] = rng.random(n) < 0.1
    df["organism"] = "sp"
    df["profile"] = "arthropod"
    df["chrom"] = np.where(np.arange(n) % 2, "c1", "c2")
    df["repro_frac"] = rng.random(n)
    return df


def _lines(out):
    return [json.loads(l) for l in (out / "folds.jsonl").read_text().splitlines()]


def test_cv_end_to_end_tiny_and_resumable(tmp_path):
    t = _tiny_table()
    n1 = cv.run_cv(t, "arthropod", ["logreg", "rf"], ["main"], 1, tmp_path, n_jobs=1)
    lines = _lines(tmp_path)
    assert n1 == len(lines) == 6                  # (2 algorithms + rules) x 1 seed x 2 units
    main = [l for l in lines if l["algorithm"] != "rules"]
    assert len(main) == 4 and all("error" not in l for l in lines)
    assert {l["unit"] for l in lines} == {"sp:c1", "sp:c2"}
    assert all(l["n_train"] + l["n_val"] == len(t) for l in lines)
    assert all(l["metrics"]["ece"] is not None for l in main)
    assert all("complete_rules_metrics" in l for l in lines)
    assert cv.run_cv(t, "arthropod", ["logreg", "rf"], ["main"], 1, tmp_path, n_jobs=1) == 0
    assert len(_lines(tmp_path)) == 6
    # a new seed only adds the missing combinations
    assert cv.run_cv(t, "arthropod", ["logreg", "rf"], ["main"], 2, tmp_path, n_jobs=1) == 4
    # and selection runs on the output without crashing
    select(tmp_path / "folds.jsonl", tmp_path / "selection.md")
    assert (tmp_path / "selection.md").read_text().startswith("# ML classifier selection")


def test_cv_failure_is_recorded_and_retried(tmp_path, monkeypatch):
    t = _tiny_table()
    real = cv.REGISTRY["logreg"]
    monkeypatch.setitem(cv.REGISTRY, "logreg",
                        lambda s, f: (_ for _ in ()).throw(RuntimeError("bad")))
    cv.run_cv(t, "arthropod", ["logreg"], ["main"], 1, tmp_path, n_jobs=1)
    errs = [l for l in _lines(tmp_path) if "error" in l]
    assert len(errs) == 2 and "bad" in errs[0]["error"] and "metrics" not in errs[0]
    with pytest.raises(SelectionError):
        select(tmp_path / "folds.jsonl", tmp_path / "s.md")
    monkeypatch.setitem(cv.REGISTRY, "logreg", real)
    assert cv.run_cv(t, "arthropod", ["logreg"], ["main"], 1, tmp_path, n_jobs=1) == 2
    select(tmp_path / "folds.jsonl", tmp_path / "s.md")      # resolved by the retry


def test_xspecies_folds():
    n = 40
    t = pd.DataFrame({"organism": ["Arabidopsis thaliana"] * n + ["Oryza sativa"] * n,
                      "chrom": ["c1", "c2"] * n, "profile": "plant"})
    units = [u for u, _, _ in cv.folds_for(t, "plant")]
    assert "xspecies:Oryza sativa" in units and "xspecies:Arabidopsis thaliana" in units
    u, tr, va = [f for f in cv.folds_for(t, "plant") if f[0] == "xspecies:Oryza sativa"][0]
    assert set(t.loc[tr, "organism"]) == {"Arabidopsis thaliana"}
    assert set(t.loc[va, "organism"]) == {"Oryza sativa"}


def test_select_class_alignment_ignores_fold_without_rules_value(tmp_path):
    rows = []
    for i, u in enumerate(UNITS):
        m = _m(0.5, amb=0.4, gap=None, ctl=None)
        if i == 0:
            m["per_class"]["siRNA"]["precision"] = None      # rules undefined on this fold
        rows.append(_line("rules", u, 0, m))
        a = _m(0.9 + 0.001 * i, prec=0.9)
        if i == 0:
            a["per_class"]["siRNA"]["precision"] = 0.0       # extreme, must be ignored
        rows.append(_line("rf", u, 0, a))
    path = tmp_path / "folds.jsonl"
    path.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    assert "**Winner: rf**" in select(path, tmp_path / "s.md")


def test_select_xspecies_excluded_from_eligibility(tmp_path):
    path = _write(tmp_path, {"rf": ([0.9 + 0.001 * i for i in range(10)], {})})
    with open(path, "a") as fh:
        fh.write(json.dumps(_line("rules", "xspecies:Oryza sativa", 0, _m(0.5, amb=0.4))) + "\n")
        fh.write(json.dumps(_line("rf", "xspecies:Oryza sativa", 0,
                                  _m(0.0, prec=0.0, rec=0.0, amb=0.99, gap=-1))) + "\n")
    text = select(path, tmp_path / "s.md")
    assert "**Winner: rf**" in text and "Cross-species folds" in text


def test_select_refuses_incomplete_seeds_or_units(tmp_path):
    path = _write(tmp_path, {"rf": ([0.8] * 10, {})})
    with pytest.raises(SelectionError, match="seed 1"):
        select(path, tmp_path / "s.md", seeds=2)
    lines = path.read_text().splitlines()
    path.write_text("\n".join(l for l in lines if '"u3"' not in l or "rf" not in l) + "\n")
    with pytest.raises(SelectionError, match="rf: seed 0 unit u3"):
        select(path, tmp_path / "s.md")


def test_read_lines_skips_truncated_final_line(tmp_path, caplog):
    p = tmp_path / "folds.jsonl"
    p.write_text(json.dumps({"a": 1}) + "\n" + '{"profile": "arth')
    assert cv.read_lines(p) == [{"a": 1}]
    p.write_text('{"bad\n' + json.dumps({"a": 1}) + "\n")
    with pytest.raises(json.JSONDecodeError):
        cv.read_lines(p)
