"""Optional ML classifier for `classify --classifier ml`.

Imported only in that mode, so the default rules mode needs neither scikit-learn nor
joblib. A model is a joblib file plus <name>.schema.json (feature list, classes, abstention
threshold, provenance); features come from mlfeatures.derive, the same code the benchmark
trained on.
"""
import json
from pathlib import Path

import mlfeatures

AMBIGUOUS = "ambiguous"
_NEEDS = "--classifier ml needs scikit-learn (pip install scikit-learn xgboost lightgbm)"


def _attr_names():
    return [*(n for n in mlfeatures.FEATURES if n != "dicer_len"), "dicer_call"]


def check_schema(schema_features):
    got, want = list(schema_features), list(mlfeatures.FEATURES)
    if got == want:
        return
    missing = [f for f in want if f not in got]
    extra = [f for f in got if f not in want]
    msg = f"model schema does not match mlfeatures.FEATURES: missing {missing}, extra {extra}"
    if not missing and not extra:
        msg += " (same names, different order)"
    raise ValueError(msg)


class MLDecider:
    def __init__(self, model, schema):
        self.model, self.schema = model, schema
        self.threshold = float(schema["threshold"])
        self.tag = f"ml:{schema['algorithm']}"

    def decide_all(self, records) -> None:
        check_schema(self.schema["features"])
        todo = [r for r in records if r.cls != "miRNA"]
        if not todo:
            return
        import numpy as np
        import pandas as pd
        names = _attr_names()
        table = pd.DataFrame([mlfeatures.derive({n: getattr(r, n, None) for n in names})
                              for r in todo], columns=list(mlfeatures.FEATURES))
        table["rule_cls"] = AMBIGUOUS
        table["carved"] = False
        proba = self.model.predict_proba(table)
        p_real = np.asarray(self.model.p_real(table), float)
        arr = proba.to_numpy()
        for i, r in enumerate(todo):
            if arr.shape[1]:
                conf = float(arr[i].max())
                cls = str(proba.columns[int(arr[i].argmax())])
            else:
                conf, cls = 0.0, AMBIGUOUS
            if conf < self.threshold:
                cls = AMBIGUOUS
            r.cls, r.confidence, r.review = cls, conf, cls == AMBIGUOUS
            r.p_real = float(min(max(p_real[i], 0.0), 1.0))
            r.evidence = self.tag
            r.score, r.class_scores = None, {}


def load(path) -> MLDecider:
    try:
        import joblib
        import sklearn  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(_NEEDS) from exc
    path = Path(path)
    schema_path = path.with_suffix(".schema.json")
    for p in (path, schema_path):
        if not p.exists():
            raise FileNotFoundError(f"--classifier ml: {p} not found")
    schema = json.loads(schema_path.read_text())
    return MLDecider(joblib.load(path), schema)
