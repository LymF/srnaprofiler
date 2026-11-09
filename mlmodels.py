"""Algorithm registry: every model is fitted inside a training fold and exposes the same API.

REGISTRY[name](seed, features) -> Model. Model.fit/predict_proba/p_real take the full table
slice (DataFrame); feature columns are selected internally (the `rules` entry needs rule_cls).

Scaling guards: real tables reach ~90k rows per profile, so hier, hdbscan, labelprop,
labelspread and svm fit on a seeded subsample of at most MAX_FIT_ROWS rows (labelled rows
always kept) and predict for all rows (1-nearest-neighbour cluster lookup for hier/hdbscan,
the estimator's own prediction otherwise). Carved miRNA rows (carved == True) are never fitted.

miRNA comes only from miRDeep2 (user decision 2026-11-07): ML_EXCLUDED_CLASSES lists the classes
the ML never outputs. Loci labelled with one of them stay PU positives (they are real) but take
the class target "other", as do non-carved control loci (`is_control`): `classes_` may contain
"other", and a decision with "other" on top is an abstention (decider.RuleScorer.choose).

p_real: supervised when the training slice has control rows (`is_control`, not carved, not
annotated): a binary model, annotated loci (any label) vs controls, fitted inside the fold with
the same family as the PU base below; `real_mode_` is "supervised". hdbscan uses the share of
annotated loci among the k nearest annotated/control training loci. Without controls the PU
score below is used (`real_mode_` = "pu"); p_real_pu always gives the PU score, for reporting.

PU score: PU classifier (positive = any label) with Elkan & Noto correction, c = mean score of a
seeded 10% held-out share of the positives, p = clip(g(x)/c, 0, 1). The PU base is logistic
regression for kmeans/gmm/hier (clustering has no supervised estimator); every supervised and
semi-supervised model uses its own estimator family as the base (logreg, rf, xgb, lgbm,
calibrated SVC, selftrain's gradient boosting; labelprop/labelspread propagate positive vs a
same-size unlabelled sample), so p_real differs between algorithms. hdbscan instead uses 1 - share of noise-labelled training loci among the k=15
nearest training loci.
"""
import numpy as np
import pandas as pd
from sklearn.calibration import CalibratedClassifierCV
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
ML_EXCLUDED_CLASSES = ("miRNA",)
OTHER = "other"   # target for excluded-class and control loci; never an output class


def _prep(features):
    return Pipeline([("imp", SimpleImputer(strategy="median", keep_empty_features=True)),
                     ("sc", StandardScaler())])


def _usable(X, y):
    keep = ~X["carved"].to_numpy(dtype=bool) if "carved" in X else np.ones(len(X), bool)
    return keep, pd.Series(y).reset_index(drop=True).where(keep, None)


STRATA_FLOOR = 50


def _subsample(n, labelled, seed, strata=None):
    """Indices (sorted) of at most ~MAX_FIT_ROWS rows; labelled rows kept when they fit.
    When labelled rows alone exceed the cap they are sampled stratified by `strata` (the
    label per row): every class keeps at least min(count, STRATA_FLOOR) rows."""
    if n <= MAX_FIT_ROWS:
        return np.arange(n)
    rng = np.random.default_rng(seed)
    lab = np.flatnonzero(labelled)
    if len(lab) > MAX_FIT_ROWS:
        return np.sort(_stratified(lab, np.asarray(strata)[lab], MAX_FIT_ROWS, rng))
    rest = np.flatnonzero(~labelled)
    room = MAX_FIT_ROWS - len(lab)
    return np.sort(np.concatenate([lab, rng.choice(rest, min(room, len(rest)), replace=False)]))


def _stratified(idx, groups, budget, rng):
    groups = np.asarray(groups)
    chosen = []
    for g in np.unique(groups):
        members = idx[groups == g]
        keep = min(len(members), STRATA_FLOOR)
        chosen.append(rng.choice(members, keep, replace=False))
    taken = np.concatenate(chosen)
    rest = np.setdiff1d(idx, taken)
    room = budget - len(taken)
    if room > 0 and len(rest):
        taken = np.concatenate([taken, rng.choice(rest, min(room, len(rest)), replace=False)])
    return taken


def _normalise(P):
    P = np.nan_to_num(np.asarray(P, float))
    s = P.sum(axis=1, keepdims=True)
    return np.where(s > 0, P / np.where(s > 0, s, 1), 1.0 / P.shape[1])


def _col1(est, Z):
    """Probability of class 1 from a fitted binary estimator (any classes_ ordering)."""
    P = est.predict_proba(Z)
    cls = list(np.asarray(est.classes_).astype(int))
    return P[:, cls.index(1)] if 1 in cls else np.zeros(len(Z))


class _PU:
    """Elkan & Noto positive-unlabelled real score.

    The binary base is the algorithm's own estimator family (ruling R20), supplied as
    `build(s) -> unfitted estimator`. mode "prop" (label propagation/spreading) fits that
    estimator on positives (1), an equally sized seeded sample of unlabelled rows (0) and the
    remaining unlabelled rows (-1). `cap` limits the fit to MAX_FIT_ROWS rows (positives kept).
    """

    def __init__(self, seed, build=None, mode="std", cap=False):
        self.seed, self.mode, self.cap = seed, mode, cap
        self.build = build or (lambda s: LogisticRegression(max_iter=1000, random_state=seed))
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
        if self.cap and len(fit_idx) > MAX_FIT_ROWS:   # stratified on positive vs unlabelled
            sub = _stratified(np.arange(len(fit_idx)), s[fit_idx], MAX_FIT_ROWS, rng)
            fit_idx = fit_idx[np.sort(sub)]
        sf = s[fit_idx]
        if self.mode == "prop":
            p_i = np.flatnonzero(sf == 1)
            u_i = np.flatnonzero(sf == 0)
            neg = rng.choice(u_i, min(len(p_i), len(u_i)), replace=False)
            y = np.full(len(fit_idx), -1)
            y[p_i], y[neg] = 1, 0
            sf = y
        self.model = self.build(sf).fit(Z[fit_idx], sf)
        if len(held):
            c = float(_col1(self.model, Z[held]).mean())
            self.c = c if c > 0 else 1.0
        return self

    def score(self, Z):
        if self.const is not None:
            return np.full(len(Z), self.const)
        return np.clip(_col1(self.model, Z) / self.c, 0.0, 1.0)


class _SupReal:
    """Supervised real score: annotated loci (1) vs negative-control loci (0).

    `build(y) -> unfitted binary estimator` is the algorithm's own family (ruling R20), the same
    builder the PU base uses; `cap` limits the fit to MAX_FIT_ROWS rows (stratified)."""

    def __init__(self, seed, build, cap=False):
        self.seed, self.build, self.cap = seed, build, cap

    def fit(self, Z, pos, neg):
        idx = np.flatnonzero(pos | neg)
        y = pos.astype(int)
        if self.cap and len(idx) > MAX_FIT_ROWS:
            idx = np.sort(_stratified(idx, y[idx], MAX_FIT_ROWS,
                                      np.random.default_rng(self.seed)))
        self.est = self.build(y[idx]).fit(Z[idx], y[idx])
        return self

    def score(self, Z):
        return _col1(self.est, Z)


class _KnnReal:
    """hdbscan's supervised real score: share of annotated loci among the K nearest
    annotated-or-control training loci (density clusterers have no estimator to supervise;
    the PU version used the noise share, which says nothing about controls)."""

    def __init__(self, seed, cap=True):
        self.seed, self.cap = seed, cap

    def fit(self, Z, pos, neg):
        idx = np.flatnonzero(pos | neg)
        if len(idx) > MAX_FIT_ROWS:
            idx = np.sort(_stratified(idx, pos[idx].astype(int), MAX_FIT_ROWS,
                                      np.random.default_rng(self.seed)))
        self.y = pos[idx].astype(float)
        self.nn = NearestNeighbors(n_neighbors=min(K_NEIGHBOURS, len(idx))).fit(Z[idx])
        return self

    def score(self, Z):
        _, ind = self.nn.kneighbors(Z)
        return self.y[ind].mean(axis=1)


class _Model:
    kind = "sup"
    real_mode_ = "pu"

    def __init__(self, name, seed, features, make):
        self.name, self.seed, self.features, self.make = name, seed, list(features), make

    def __getstate__(self):
        # `make` is a lambda/closure and cannot be pickled; a loaded model only predicts.
        return {**self.__dict__, "make": None}

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.loaded_ = True

    # --- API
    def fit(self, X, y):
        if getattr(self, "loaded_", False):
            raise RuntimeError("this model was loaded from an artifact and cannot be refitted; "
                               "build a new one with REGISTRY[name](seed, features)")
        keep, y = _usable(X, y)
        self.prep_ = _prep(self.features)
        Z = self.prep_.fit_transform(X[self.features].to_numpy(float)[keep])
        y = y[keep].reset_index(drop=True)
        real = y.notna().to_numpy()                      # PU positives: any label
        # class targets: excluded classes (miRNA) and non-carved control loci become the
        # explicit "other" target, so the model can abstain instead of forcing a class
        ctrl = (X["is_control"].to_numpy(dtype=bool)[keep] if "is_control" in X
                else np.zeros(len(y), bool))
        other = y.isin(ML_EXCLUDED_CLASSES).to_numpy() | (ctrl & ~real)
        y = y.where(~other, OTHER)
        lab = y.notna().to_numpy()
        self.classes_ = sorted(y.dropna().unique())
        self.pu_ = self._pu().fit(Z, real)
        neg = ctrl & ~real                               # negative controls (non-carved)
        self.real_mode_ = "pu"
        self.sup_ = None
        if neg.any() and real.any():                     # supervised real needs both sides
            self.sup_ = self._real_sup().fit(Z, real, neg)
            self.real_mode_ = "supervised"
        if self.classes_:
            self._fit_model(Z, y, lab)
        return self

    def _pu(self):
        """PU real-score base; clustering models have no supervised estimator, so logistic."""
        return _PU(self.seed)

    def _real_sup(self):
        """Supervised real-score base: the family's own estimator (logistic for clustering)."""
        return _SupReal(self.seed, lambda y: LogisticRegression(
            max_iter=1000, random_state=self.seed))

    def predict_proba(self, X):
        Z = self.prep_.transform(X[self.features].to_numpy(float))
        if not self.classes_:
            return pd.DataFrame(index=X.index, columns=[], dtype=float)
        return pd.DataFrame(self._proba(Z), index=X.index, columns=self.classes_)

    def p_real(self, X):
        """Supervised score (annotated vs control) when the training slice had controls,
        else the PU score (`real_mode_` says which)."""
        if self.sup_ is not None:
            return self.sup_.score(self.prep_.transform(X[self.features].to_numpy(float)))
        return self.p_real_pu(X)

    def p_real_pu(self, X):
        """Elkan & Noto PU score, kept for reporting."""
        return self.pu_.score(self.prep_.transform(X[self.features].to_numpy(float)))

    def _enc(self, y):
        idx = {c: i for i, c in enumerate(self.classes_)}
        return y.map(idx).to_numpy()


class _Supervised(_Model):
    cap = False

    def _fit_model(self, Z, y, lab):
        idx = np.flatnonzero(lab)
        if self.cap and len(idx) > MAX_FIT_ROWS:
            idx = np.sort(_stratified(idx, y.to_numpy()[idx], MAX_FIT_ROWS,
                                      np.random.default_rng(self.seed)))
        self.fit_rows_ = idx
        yy = self._enc(y)[idx].astype(int)
        self.const_ = None
        if len(self.classes_) < 2:
            self.const_ = 1.0
            return
        self.est_ = self._build(yy).fit(Z[idx], yy)

    def _build(self, yy):
        return self.make(self.seed)

    def _pu(self):
        return _PU(self.seed, self._build, "std", self.cap)

    def _real_sup(self):
        return _SupReal(self.seed, self._build, self.cap)

    def _proba(self, Z):
        if self.const_ is not None:
            return np.ones((len(Z), 1))
        P = np.zeros((len(Z), len(self.classes_)))
        P[:, np.asarray(self.est_.classes_).astype(int)] = self.est_.predict_proba(Z)
        return _normalise(P)


class _OneHotSVC:
    """Fallback when a class has < 2 rows (too few to calibrate): hard SVC predictions."""

    def __init__(self, seed):
        self.svc = SVC(random_state=seed)

    def fit(self, Z, y):
        self.svc.fit(Z, y)
        self.classes_ = self.svc.classes_
        return self

    def predict_proba(self, Z):
        return (self.svc.predict(Z)[:, None] == self.classes_[None, :]).astype(float)


class _SVM(_Supervised):
    """SVC calibrated with CalibratedClassifierCV (sigmoid, 3-fold, single refit-free ensemble);
    fitted on at most MAX_FIT_ROWS labelled rows, stratified by class."""
    cap = True

    def _build(self, yy):
        cv = min(3, int(np.bincount(yy).min()) if len(yy) else 0)
        if cv < 2:
            return _OneHotSVC(self.seed)
        return CalibratedClassifierCV(SVC(random_state=self.seed), method="sigmoid", cv=cv,
                                      ensemble=False)


class _SemiSup(_Model):
    """labelprop/labelspread/selftrain: labelled rows kept, unlabelled subsampled to the cap."""

    def _pu(self):
        if self.name in ("labelprop", "labelspread"):
            return _PU(self.seed, lambda s: self.make(self.seed), "prop", True)
        return _PU(self.seed, lambda s: _hgb(self.seed), "std", False)   # selftrain base

    def _real_sup(self):
        if self.name in ("labelprop", "labelspread"):
            return _SupReal(self.seed, lambda y: self.make(self.seed), True)
        return _SupReal(self.seed, lambda y: _hgb(self.seed), True)

    def _fit_model(self, Z, y, lab):
        idx = _subsample(len(Z), lab, self.seed, y.to_numpy())
        self.fit_rows_ = idx
        yy = np.where(lab, self._enc(y), -1).astype(int)[idx]
        self.const_ = None
        if len(set(yy[yy >= 0])) < 2:
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
        idx = _subsample(len(Z), lab, self.seed, y.to_numpy())
        self.fit_rows_ = idx
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

    def _real_sup(self):
        return _KnnReal(self.seed)

    def p_real_pu(self, X):
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
    "svm": _factory(_SVM, None, "svm"),
    "rules": lambda seed, features: _Rules(seed, features),
}
