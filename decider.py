"""Pluggable decider. v1: RuleScorer (transparent scoring model)."""
import math
import yaml


def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def _sigmoid(x, x0=1.65, scale=1.0):
    return 1.0 / (1.0 + math.exp(-(x - x0) / scale))


def _window(size, lo, hi, edge=3.0):
    if size is None:
        return 0.0
    if lo <= size <= hi:
        return 1.0
    d = lo - size if size < lo else size - hi
    return max(0.0, 1.0 - d / edge)


def _mid(frac_top, sigma=0.15):
    if frac_top is None:
        return 0.0
    return math.exp(-((frac_top - 0.5) ** 2) / (2 * sigma ** 2))


class RuleScorer:
    def __init__(self, config):
        self.cfg = config
        # Configurable shape parameters (defaults = previous behavior).
        self.shape = (config.get("decision", {}) or {}).get("shape", {}) or {}

    def _size(self, record):
        return int(record.dicer_call) if record.dicer_call != "N" else None

    def _pp_signal(self, record):
        # Ping-pong strength from the bounded 10-nt overlap fraction (comparable
        # across loci, unlike the raw z). Linear ramp: 0 at frac 0, saturating to 1
        # at `pingpong_frac_full` (default 0.1 = 10% of opposite-strand pairs overlap by 10).
        full = self.shape.get("pingpong_frac_full", 0.1)
        return min(record.pingpong_frac / full, 1.0) if full > 0 else 0.0

    def _score_class(self, name, spec, record):
        size = self._size(record)
        edge = self.shape.get("size_edge", 3.0)
        sigma = self.shape.get("fractop_sigma", 0.15)
        # Size gate for the ping-pong terms: a 10-nt 5'-5' overlap only carries
        # ping-pong evidence where the size is already compatible with the class.
        # Dual-strand 21-nt loci (canonical endo-siRNA) hit that overlap by geometry
        # alone; without the gate, w_pingpong (the heaviest weight) plus the siRNA
        # penalty flipped them to piRNA against an unambiguous size call.
        gate = (_window(size, spec["size"]["lo"], spec["size"]["hi"], edge=edge)
                if "size" in spec else 1.0)
        s = 0.0
        terms = {}
        if "w_size" in spec:
            t = spec["w_size"] * _window(size, spec["size"]["lo"], spec["size"]["hi"], edge=edge)
            s += t; terms["size"] = round(t, 3)
        if "w_pingpong" in spec:
            t = spec["w_pingpong"] * self._pp_signal(record) * gate
            s += t; terms["pingpong"] = round(t, 3)
        if "w_1u" in spec:
            t = spec["w_1u"] * record.pct_1u
            s += t; terms["1u"] = round(t, 3)
        if "w_fractop_mid" in spec:
            t = spec["w_fractop_mid"] * _mid(record.frac_top, sigma=sigma)
            s += t; terms["fractop_mid"] = round(t, 3)
        if "w_pingpong_penalty" in spec:
            t = spec["w_pingpong_penalty"] * self._pp_signal(record) * gate
            s -= t; terms["pingpong_penalty"] = round(-t, 3)
        if "w_phase" in spec:
            # phase_score is unbounded (ratio*number*log(abun)); saturate to [0,1].
            t = spec["w_phase"] * _sigmoid(
                record.phase_score,
                x0=spec.get("phase_x0", 5.0), scale=spec.get("phase_scale", 3.0))
            s += t; terms["phase"] = round(t, 3)
        return s, terms

    def decide(self, record):
        scores = {}
        evidence = {}
        for name, spec in self.cfg["classes"].items():
            scores[name], evidence[name] = self._score_class(name, spec, record)
        record.class_scores = {k: round(v, 3) for k, v in scores.items()}
        # softmax for confidence
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
