"""Decisor plugável. v1: RuleScorer (modelo de score transparente)."""
import math
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def _sigmoid(x, x0=1.65, scale=1.0):
    return 1.0 / (1.0 + math.exp(-(x - x0) / scale))


def _window(size, lo, hi):
    if size is None:
        return 0.0
    if lo <= size <= hi:
        return 1.0
    d = lo - size if size < lo else size - hi
    return max(0.0, 1.0 - d / 3.0)


def _mid(frac_top):
    if frac_top is None:
        return 0.0
    return math.exp(-((frac_top - 0.5) ** 2) / (2 * 0.15 ** 2))


class RuleScorer:
    def __init__(self, config):
        self.cfg = config

    def _size(self, record):
        return int(record.dicer_call) if record.dicer_call != "N" else None

    def _score_class(self, name, spec, record):
        size = self._size(record)
        s = 0.0
        terms = {}
        if "w_size" in spec:
            t = spec["w_size"] * _window(size, spec["size"]["lo"], spec["size"]["hi"])
            s += t; terms["size"] = round(t, 3)
        if "w_pingpong" in spec:
            t = spec["w_pingpong"] * _sigmoid(record.pingpong_z)
            s += t; terms["pingpong"] = round(t, 3)
        if "w_1u" in spec:
            t = spec["w_1u"] * record.pct_1u
            s += t; terms["1u"] = round(t, 3)
        if "w_fractop_mid" in spec:
            t = spec["w_fractop_mid"] * _mid(record.frac_top)
            s += t; terms["fractop_mid"] = round(t, 3)
        if "w_pingpong_penalty" in spec:
            t = spec["w_pingpong_penalty"] * _sigmoid(record.pingpong_z)
            s -= t; terms["pingpong_penalty"] = round(-t, 3)
        return s, terms

    def decide(self, record):
        scores = {}
        evidence = {}
        for name, spec in self.cfg["classes"].items():
            scores[name], evidence[name] = self._score_class(name, spec, record)
        record.class_scores = {k: round(v, 3) for k, v in scores.items()}
        # softmax para confiança
        ordered = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
        top_name, top_score = ordered[0]
        exps = {k: math.exp(v) for k, v in scores.items()}
        z = sum(exps.values())
        conf = exps[top_name] / z if z else 0.0
        margin = top_score - (ordered[1][1] if len(ordered) > 1 else 0.0)
        dec = self.cfg["decision"]
        record.cls = top_name
        record.confidence = round(conf, 3)
        record.evidence = "; ".join(f"{k}={v}" for k, v in evidence[top_name].items())
        record.review = (
            conf < dec["review_confidence"]
            or margin < dec["min_margin"]
            or record.reads < dec["min_reads"]
            or record.dicer_call == "N"
        )
        return record
