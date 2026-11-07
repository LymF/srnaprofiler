"""Algorithm registry: every model is fitted inside a training fold and exposes the same API.

REGISTRY[name](seed, features) -> Model. Model.fit/predict_proba/p_real take the full table
slice (DataFrame); feature columns are selected internally (the `rules` entry needs rule_cls).

Scaling guards: real tables reach ~90k rows per profile, so hier, hdbscan, labelprop,
labelspread and svm fit on a seeded subsample of at most MAX_FIT_ROWS rows (labelled rows
always kept) and predict for all rows (1-nearest-neighbour cluster lookup for hier/hdbscan,
the estimator's own prediction otherwise). Carved miRNA rows (carved == True) are never fitted.

p_real: PU classifier (positive = any label) with Elkan & Noto correction, c = mean score of a
seeded 10% held-out share of the positives, p = clip(g(x)/c, 0, 1). The PU base is logistic
regression. hdbscan instead uses 1 - share of noise-labelled training loci among the k=15
nearest training loci.
"""
import numpy as np
import pandas as pd
from sklearn.cluster import HDBSCAN, AgglomerativeClustering, KMeans
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.semi_supervised import LabelPropagation, LabelSpreading, SelfTrainingClassifier
from sklearn.svm import SVC

MAX_FIT_ROWS = 20000
K_NEIGHBOURS = 15
PU_HOLDOUT = 0.1


def _prep(features):
    return Pipeline([("imp", SimpleImputer(strategy="median", keep_empty_features=True)),
                     ("sc", StandardScaler())])


def _usable(X, y):
    keep = ~X["carved"].to_numpy(dtype=bool) if "carved" in X else np.ones(len(X), bool)
    return keep, pd.Series(y).reset_index(drop=True).where(keep, None)


def _subsample(n, labelled, seed):
    """Indices (sorted) of at most MAX_FIT_ROWS rows; labelled rows always kept."""
    if n <= MAX_FIT_ROWS:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    lab = np.flatnonzero(labelled)
    rest = np.flatnonzero(~labelled)
    room = max(MAX_FIT_ROWS - len(lab), 0)
    return np.sort(np.concatenate([lab, rng.choice(rest, min(room, len(rest)), replace=False)]))


def _normalise(P):
    P = np.nan_to_num(np.asarray(P, float))
    s = P.sum(axis=1, keepdims=True)
    return np.where(s > 0, P / np.where(s > 0, s, 1), 1.0 / P.shape[1])


class _PU:
    """Elkan & Noto positive-unlabelled real score."""

    def __init__(self, seed):
        self.seed = seed
        self.model = None
        self.c = 1.0
        self.const = None

    def fit(self, Z, positive):
        pos = np.flatnonzero(positive)
        if len(pos) == 0:
            self.const = 0.0
            return self
        if len(pos) == len(positive):
            self.const = 1.0
            return self
        rng = np.random.default_rng(self.seed)
        held = np.array([], int)
        if len(pos) >= 2:
            held = rng.choice(pos, max(1, int(round(PU_HOLDOUT * len(pos)))), replace=False)
        s = positive.astype(int).copy()
        s[held] = 0
        fit_idx = np.setdiff1d(np.arange(len(Z)), held)
        self.model = LogisticRegression(max_iter=1000, random_state=self.seed)
        self.model.fit(Z[fit_idx], s[fit_idx])
        if len(held):
            c = float(self.model.predict_proba(Z[held])[:, 1].mean())
            self.c = c if c > 0 else 1.0
        return self

    def score(self, Z):
        if self.const is not None:
            return np.full(len(Z), self.const)
        return np.clip(self.model.predict_proba(Z)[:, 1] / self.c, 0.0, 1.0)


class _Model:
    kind = "sup"

    def __init__(self, name, seed, features, make):
        self.name, self.seed, self.features, self.make = name, seed, list(features), make

    # --- API
    def fit(self, X, y):
        keep, y = _usable(X, y)
        self.prep_ = _prep(self.features)
        Z = self.prep_.fit_transform(X[self.features].to_numpy(float)[keep])
        y = y[keep].reset_index(drop=True)
        lab = y.notna().to_numpy()
        self.classes_ = sorted(y.dropna().unique())
        self.pu_ = _PU(self.seed).fit(Z, lab)
        if self.classes_:
            self._fit_model(Z, y, lab)
        return self

    def predict_proba(self, X):
        Z = self.prep_.transform(X[self.features].to_numpy(float))
        if not self.classes_:
            return pd.DataFrame(index=X.index, columns=[], dtype=float)
        return pd.DataFrame(self._proba(Z), index=X.index, columns=self.classes_)

    def p_real(self, X):
        return self.pu_.score(self.prep_.transform(X[self.features].to_numpy(float)))

    def _enc(self, y):
        idx = {c: i for i, c in enumerate(self.classes_)}
        return y.map(idx).to_numpy()


class _Supervised(_Model):
    cap = False

    def _fit_model(self, Z, y, lab):
        idx = np.flatnonzero(lab)
        if self.cap and len(idx) > MAX_FIT_ROWS:
            idx = np.sort(np.random.default_rng(self.seed).choice(idx, MAX_FIT_ROWS, replace=False))
        yy = self._enc(y)[idx].astype(int)
        self.const_ = None
        if len(self.classes_) < 2:
            self.const_ = 1.0
            return
        self.est_ = self.make(self.seed).fit(Z[idx], yy)

    def _proba(self, Z):
        if self.const_ is not None:
            return np.ones((len(Z), 1))
        return _normalise(self.est_.predict_proba(Z))


class _SVM(_Supervised):
    cap = True


class _SemiSup(_Model):
    """labelprop/labelspread/selftrain: labelled rows kept, unlabelled subsampled to the cap."""

    def _fit_model(self, Z, y, lab):
        idx = _subsample(len(Z), lab, self.seed)
        yy = np.where(lab, self._enc(y), -1).astype(int)[idx]
        self.const_ = None
        if len(self.classes_) < 2:
            self.const_ = 1.0
            return
        self.est_ = self.make(self.seed).fit(Z[idx], yy)

    def _proba(self, Z):
        if self.const_ is not None:
            return np.ones((len(Z), 1))
        P = np.zeros((len(Z), len(self.classes_)))
        P[:, self.est_.classes_.astype(int)] = self.est_.predict_proba(Z)
        return _normalise(P)


class _Unsup(_Model):
    kind = "unsup"

    def _fit_model(self, Z, y, lab):
        idx = _subsample(len(Z), lab, self.seed)
        self.Zfit_ = Z[idx]
        k = len(self.classes_) + 2
        self.cl_ = self._cluster(self.Zfit_, k)
        self.nn_ = NearestNeighbors(n_neighbors=min(K_NEIGHBOURS, len(idx))).fit(self.Zfit_)
        lab_i = lab[idx]
        yy = y.to_numpy()[idx]
        self.map_ = {}
        for c in np.unique(self.cl_):
            if c < 0:
                continue
            votes = pd.Series(yy[(self.cl_ == c) & lab_i])
            if len(votes):
                self.map_[c] = votes.value_counts().sort_index().idxmax()

    def _cluster(self, Z, k):
        raise NotImplementedError

    def _assign(self, Z):
        """Cluster id per row; default 1-nearest-neighbour among the fitted points."""
        _, ind = self.nn_.kneighbors(Z, n_neighbors=1)
        return self.cl_[ind[:, 0]]

    def _proba(self, Z):
        P = np.zeros((len(Z), len(self.classes_)))
        col = {c: i for i, c in enumerate(self.classes_)}
        for r, c in enumerate(self._assign(Z)):
            cls = self.map_.get(c)
            if cls is None:
                P[r] = 1.0 / len(self.classes_)
            else:
                P[r, col[cls]] = 1.0
        return P


class _KMeans(_Unsup):
    def _cluster(self, Z, k):
        self.km_ = KMeans(n_clusters=min(k, len(Z)), n_init=5, random_state=self.seed).fit(Z)
        return self.km_.labels_

    def _assign(self, Z):
        return self.km_.predict(Z)


class _GMM(_Unsup):
    def _cluster(self, Z, k):
        self.gm_ = GaussianMixture(n_components=min(k, len(Z)), reg_covar=1e-3,
                                   random_state=self.seed).fit(Z)
        return self.gm_.predict(Z)

    def _assign(self, Z):
        return self.gm_.predict(Z)


class _Hier(_Unsup):
    def _cluster(self, Z, k):
        return AgglomerativeClustering(n_clusters=min(k, len(Z))).fit_predict(Z)


class _HDBSCAN(_Unsup):
    def _cluster(self, Z, k):
        return HDBSCAN(min_cluster_size=max(5, min(15, len(Z) // 10)),
                       copy=True).fit_predict(Z)

    def p_real(self, X):
        Z = self.prep_.transform(X[self.features].to_numpy(float))
        if not self.classes_:
            return np.zeros(len(Z))
        _, ind = self.nn_.kneighbors(Z)
        return 1.0 - (self.cl_[ind] < 0).mean(axis=1)


class _Rules:
    """Baseline passthrough of rule_cls; no fitting."""

    def __init__(self, seed, features):
        pass

    def fit(self, X, y):
        keep, y = _usable(X, y)
        self.classes_ = sorted(y.dropna().unique())
        return self

    def predict_proba(self, X):
        P = pd.DataFrame(0.0, index=X.index, columns=self.classes_)
        rc = X["rule_cls"].to_numpy(dtype=object)
        for cls in self.classes_:
            P[cls] = (rc == cls).astype(float)
        none = P.sum(axis=1) == 0
        P.loc[none, :] = 1.0 / max(len(self.classes_), 1)
        return P

    def p_real(self, X):
        return (X["rule_cls"].to_numpy(dtype=object) != "ambiguous").astype(float)


def _lazy(module, attr, **kw):
    def make(seed):
        try:
            mod = __import__(module, fromlist=[attr])
        except Exception as exc:  # e.g. missing libomp raises XGBoostError/OSError
            raise ImportError(f"{module} unavailable: {exc}") from exc
        return getattr(mod, attr)(random_state=seed, **kw)
    return make


def _factory(cls, make=None, name=None, probe=None):
    def build(seed, features):
        if probe:
            probe(seed)
        return cls(name, seed, features, make)
    return build


_xgb = _lazy("xgboost", "XGBClassifier", n_estimators=100, n_jobs=1, verbosity=0)
_lgbm = _lazy("lightgbm", "LGBMClassifier", n_estimators=100, n_jobs=1, verbose=-1)


def _hgb(seed):
    return HistGradientBoostingClassifier(max_iter=50, random_state=seed)


REGISTRY = {
    "kmeans": _factory(_KMeans, name="kmeans"),
    "gmm": _factory(_GMM, name="gmm"),
    "hdbscan": _factory(_HDBSCAN, name="hdbscan"),
    "hier": _factory(_Hier, name="hier"),
    "labelprop": _factory(_SemiSup, lambda s: LabelPropagation(kernel="knn", n_neighbors=10,
                                                               max_iter=200), "labelprop"),
    "labelspread": _factory(_SemiSup, lambda s: LabelSpreading(kernel="knn", n_neighbors=10,
                                                               max_iter=200), "labelspread"),
    "selftrain": _factory(_SemiSup, lambda s: SelfTrainingClassifier(_hgb(s)), "selftrain"),
    "logreg": _factory(_Supervised, lambda s: LogisticRegression(max_iter=1000, random_state=s),
                       "logreg"),
    "rf": _factory(_Supervised, lambda s: RandomForestClassifier(
        n_estimators=100, n_jobs=1, random_state=s), "rf"),
    "xgb": _factory(_Supervised, _xgb, "xgb", probe=_xgb),
    "lgbm": _factory(_Supervised, _lgbm, "lgbm", probe=_lgbm),
    "svm": _factory(_SVM, lambda s: SVC(probability=True, random_state=s), "svm"),
    "rules": lambda seed, features: _Rules(seed, features),
}
