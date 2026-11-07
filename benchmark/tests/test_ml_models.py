import numpy as np
import pandas as pd
import pytest

pytest.importorskip("sklearn")

import mlmodels as models
from benchmark.ml import folds
from benchmark.ml.models import REGISTRY
from mlfeatures import FEATURES

CLASSES = ["miRNA", "siRNA", "piRNA"]
FEATS = FEATURES[:6] + ["ctx_te", "ctx_gene"] + ["repro_frac"]


def make_table(n=150, seed=0):
    rng = np.random.default_rng(seed)
    cls = rng.integers(0, 3, n)
    X = rng.normal(0, 0.5, (n, 6)) + cls[:, None] * 4.0
    df = pd.DataFrame(X, columns=FEATS[:6])
    df["ctx_te"] = np.nan
    df["ctx_gene"] = np.nan
    df["repro_frac"] = np.nan
    lab = np.array(CLASSES, dtype=object)[cls]
    keep = rng.random(n) < 0.3
    df["label"] = np.where(keep, lab, None)
    df["true"] = lab
    df["rule_cls"] = np.where(rng.random(n) < 0.2, "ambiguous", lab)
    df["carved"] = False
    df["organism"] = "sp"
    df["profile"] = "arthropod"
    df["chrom"] = rng.choice(["c1", "c2"], n)
    return df


def _skip_missing(name, seed=0):
    try:
        return REGISTRY[name](seed, FEATS)
    except ImportError:
        pytest.skip(f"{name}: library unavailable")


@pytest.mark.parametrize("name", sorted(REGISTRY))
def test_every_model_fits_with_nan_columns(name):
    df = make_table()
    m = _skip_missing(name)
    m.fit(df, df["label"])
    P = m.predict_proba(df)
    assert list(P.columns) == sorted(df["label"].dropna().unique())
    assert np.allclose(P.sum(axis=1), 1.0, atol=1e-6)
    p = m.p_real(df)
    assert p.shape == (len(df),) and ((p >= 0) & (p <= 1)).all()


def test_registry_has_all_names():
    assert set(REGISTRY) == {"kmeans", "gmm", "hdbscan", "hier", "labelprop", "labelspread",
                             "selftrain", "logreg", "rf", "xgb", "lgbm", "svm", "rules"}
    assert models.MAX_FIT_ROWS == 20000


def test_supervised_beats_chance_on_blobs():
    df = make_table(300)
    lab = df["label"].notna().to_numpy()
    idx = np.flatnonzero(lab)
    tr, te = idx[::2], idx[1::2]
    m = REGISTRY["rf"](0, FEATS)
    m.fit(df.iloc[tr], df["label"].iloc[tr])
    pred = m.predict_proba(df.iloc[te]).idxmax(axis=1)
    assert (pred.to_numpy() == df["true"].iloc[te].to_numpy()).mean() > 0.9


def test_carved_rows_excluded_from_fit():
    df = make_table()
    df["carved"] = False
    df.loc[df["label"] == "piRNA", "carved"] = True
    m = REGISTRY["rf"](0, FEATS).fit(df, df["label"])
    assert "piRNA" not in m.predict_proba(df).columns


def test_mapping_fitted_on_train_only():
    df = make_table()
    y = df["label"].copy()
    train = df.iloc[:100]
    y_tr = y.iloc[:100].copy()
    y_tr[y_tr == "piRNA"] = None
    for name in ("kmeans", "gmm", "hdbscan", "hier"):
        m = REGISTRY[name](0, FEATS).fit(train, y_tr)
        P = m.predict_proba(df)
        assert "piRNA" not in P.columns
        assert set(P.columns) == {"miRNA", "siRNA"}


def test_rules_passthrough():
    df = make_table()
    m = REGISTRY["rules"](0, FEATS).fit(df, df["label"])
    P = m.predict_proba(df)
    amb = (df["rule_cls"] == "ambiguous").to_numpy()
    assert np.allclose(P.sum(axis=1), 1.0)
    assert np.allclose(P[amb], 1 / P.shape[1])
    assert (P[~amb].idxmax(axis=1).to_numpy() == df["rule_cls"][~amb].to_numpy()).all()
    assert (m.p_real(df) == (~amb).astype(float)).all()


def test_chromosome_units_pool_small_contigs():
    chrom = ["c1"] * 60 + ["c2"] * 38 + ["c3"] * 2
    df = pd.DataFrame({"organism": "sp", "profile": "p", "chrom": chrom})
    u = folds.chromosome_units(df, 0.05)
    assert set(u) == {"sp:c1", "sp:c2", "sp:other"}


def test_units_unique_across_organisms():
    df = pd.DataFrame({"organism": ["a"] * 10 + ["b"] * 10, "chrom": ["1"] * 20})
    assert set(folds.chromosome_units(df)) == {"a:1", "b:1"}


def test_folds_never_mix_units():
    chrom = ["c1"] * 60 + ["c2"] * 38 + ["c3"] * 2
    df = pd.DataFrame({"organism": "sp", "profile": "p", "chrom": chrom})
    df.index = range(100, 200)
    units = folds.chromosome_units(df)
    got = list(folds.split(df))
    assert len(got) == 3
    for tr, va, unit in got:
        assert set(units[tr]).isdisjoint(set(units[va]))
        assert set(units[va]) == {unit}
        assert len(tr) + len(va) == 100


def _cap_table(n_lab_per_class, n_unlab, rare=None, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for k, c in enumerate(CLASSES):
        m = rare if (rare is not None and c == "piRNA") else n_lab_per_class
        rows += [(k, c)] * m
    rows += [(rng.integers(0, 3), None)] * n_unlab
    cls = np.array([r[0] for r in rows])
    df = pd.DataFrame(rng.normal(0, 0.5, (len(rows), 6)) + cls[:, None] * 4.0, columns=FEATS[:6])
    for c in ("ctx_te", "ctx_gene", "repro_frac"):
        df[c] = np.nan
    df["label"] = pd.Series([r[1] for r in rows], dtype=object)
    df["rule_cls"] = "siRNA"
    df["carved"] = False
    return df.sample(frac=1, random_state=seed).reset_index(drop=True)


@pytest.mark.parametrize("name", ["hier", "hdbscan", "labelprop", "svm"])
def test_cap_path_keeps_labelled_and_predicts_all(name, monkeypatch):
    monkeypatch.setattr(models, "MAX_FIT_ROWS", 30)
    df = _cap_table(4, 88)          # 12 labelled rows fit under the cap of 30
    m = REGISTRY[name](0, FEATS).fit(df, df["label"])
    rows = m.fit_rows_
    assert len(rows) <= 30
    assert set(np.flatnonzero(df["label"].notna())) <= set(rows)
    P = m.predict_proba(df)
    assert len(P) == len(df) and np.allclose(P.sum(axis=1), 1.0, atol=1e-6)


@pytest.mark.parametrize("name", ["svm", "labelprop", "hier"])
def test_cap_stratified_keeps_rare_class(name, monkeypatch):
    monkeypatch.setattr(models, "MAX_FIT_ROWS", 120)
    df = _cap_table(100, 20, rare=3)    # 203 labelled > cap, piRNA has 3 rows
    m = REGISTRY[name](0, FEATS).fit(df, df["label"])
    assert len(m.fit_rows_) <= 120 + 3
    P = m.predict_proba(df)
    assert list(P.columns) == sorted(CLASSES) and len(P) == len(df)
    assert np.allclose(P.sum(axis=1), 1.0, atol=1e-6)
    assert df["label"].iloc[m.fit_rows_].eq("piRNA").sum() == 3
