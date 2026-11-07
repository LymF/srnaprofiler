import numpy as np
import pandas as pd
import pytest

from benchmark.ml import evaluate as ev


class Fixed:
    def __init__(self, P, classes):
        self.P, self.classes = np.asarray(P, float), classes

    def predict_proba(self, table):
        return pd.DataFrame(self.P, columns=self.classes, index=table.index)


def test_macro_f1_skips_classes_without_truth():
    val = pd.DataFrame({"label": ["piRNA", "piRNA", None], "repro_frac": [np.nan] * 3})
    pred = ["piRNA", "ambiguous", "siRNA"]
    m = ev.fold_metrics(val, pred, [0.9, 0.9, 0.9], ["piRNA"] * 3)
    assert m["per_class"]["siRNA"]["precision"] is None
    assert m["per_class"]["piRNA"]["precision"] == 1.0
    assert m["macro_f1"] == pytest.approx(2 * 1.0 * 0.5 / 1.5)
    assert m["ambiguous_share"] == pytest.approx(1 / 3)
    assert m["real_repro_gap"] is None


def test_decide_abstains_below_threshold():
    t = pd.DataFrame({"carved": [False, False]})
    m = Fixed([[0.5, 0.5], [0.9, 0.1]], ["a", "b"])
    assert list(ev.decide(m, t, [], 0.6)) == ["ambiguous", "a"]


def test_carved_stays_mirna():
    t = pd.DataFrame({"carved": [True]})
    assert list(ev.decide(Fixed([[0.1, 0.9]], ["a", "b"]), t, [], 0.0)) == ["miRNA"]


def test_complete_rules_keeps_rule_calls():
    out = ev.complete_rules(["piRNA", "ambiguous"], ["siRNA", "piRNA"])
    assert list(out) == ["piRNA", "piRNA"]


def test_control_rates():
    val = pd.DataFrame({"label": [None, None], "is_control": [True, True],
                        "repro_frac": [0.8, 0.2]})
    m = ev.fold_metrics(val, ["a", "a"], [0.9, 0.1], ["piRNA", "ambiguous"])
    assert m["control_real_rate"] == 0.5
    assert m["rules_control_rate"] == 0.5
    assert m["real_repro_gap"] == pytest.approx(0.6)


def _table(n=90):
    rng = np.random.default_rng(0)
    y = np.array(["siRNA", "piRNA"])[rng.integers(0, 2, n)]
    return pd.DataFrame({"label": y, "chrom": np.repeat(["c1", "c2", "c3"], n // 3),
                         "organism": "sp", "carved": False, "rule_cls": "ambiguous"}), y


class Oracle:
    def __init__(self, perfect, seed=0):
        self.perfect, self.rng = perfect, np.random.default_rng(seed)

    def fit(self, table):
        return self

    def predict_proba(self, table):
        if self.perfect:
            P = np.stack([(table["label"] == "piRNA").to_numpy(),
                          (table["label"] == "siRNA").to_numpy()], 1).astype(float)
        else:
            P = np.tile([0.5, 0.5], (len(table), 1))
        return pd.DataFrame(P, columns=["piRNA", "siRNA"], index=table.index)


def test_choose_threshold_perfect_model_returns_zero():
    t, y = _table()
    t["rule_cls"] = np.where(np.arange(len(t)) % 2, y, "ambiguous")
    t.loc[::4, "rule_cls"] = np.where(y[::4] == "siRNA", "piRNA", "siRNA")
    assert ev.choose_threshold(lambda: Oracle(True), t, [], 0) == 0.0


def test_choose_threshold_random_model_perfect_rules_returns_fallback():
    t, y = _table()
    t["rule_cls"] = y
    assert ev.choose_threshold(lambda: Oracle(False), t, [], 0) == 0.95


def test_decide_never_mirna_on_noncarved():
    from benchmark.ml.models import REGISTRY
    from benchmark.tests.test_ml_models import FEATS, make_table, with_mirna
    df = with_mirna(make_table())
    df.loc[:9, "carved"] = True
    m = REGISTRY["logreg"](0, FEATS)
    m.fit(df, df["label"])
    out = ev.decide(m, df, FEATS, 0.0)
    assert (out[df["carved"]] == "miRNA").all()
    assert not (out[~df["carved"]] == "miRNA").any()


def test_choose_threshold_ignores_mirna_labelled_rows():
    t, y = _table()
    t["rule_cls"] = y
    base = ev.choose_threshold(lambda: Oracle(True), t, [], 0)
    assert base == 0.0
    extra = t.iloc[:30].copy()
    extra["label"] = "miRNA"
    extra["rule_cls"] = "ambiguous"
    both = pd.concat([t, extra], ignore_index=True)
    # miRNA rows get zero probability everywhere; counting them as errors would raise t
    assert ev.choose_threshold(lambda: Oracle(True), both, [], 0) == 0.0
