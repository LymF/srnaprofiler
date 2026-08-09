"""Pluggable decider. v1: RuleScorer (transparent scoring model)."""
import math

import yaml

from evidence import effect_strength
from signals import PHASE_WINDOW_CYCLES


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
    # Only names with an e_* attribute on ClusterRecord; others silently score 0.0.
    KNOWN_EVIDENCE = ("1u", "5a", "fractop_bid", "pingpong", "phase", "size")

    PHASE_MODES = ("score", "cycles")

    def __init__(self, config):
        self.cfg = config
        # Configurable shape parameters (defaults = previous behavior).
        self.shape = (config.get("decision", {}) or {}).get("shape", {}) or {}
        for name, spec in (config.get("classes") or {}).items():
            mode = spec.get("phase_on", "score")
            if mode not in self.PHASE_MODES:
                raise ValueError(
                    f"{name}: unknown phase_on {mode!r}; "
                    f"known: {', '.join(self.PHASE_MODES)}")
            for entry in spec.get("defining_evidence", []):
                for ev in (entry if isinstance(entry, list) else [entry]):
                    if ev not in self.KNOWN_EVIDENCE:
                        raise KeyError(
                            f"{name}: unknown evidence {ev!r} in "
                            f"defining_evidence; known: {', '.join(self.KNOWN_EVIDENCE)}")

    def size_strength(self, spec, record):
        """How much of the locus sits in this class's size window, as a strength.

        The only class-dependent evidence: the other e_* live on the record
        because they mean the same thing for every class, but a size only means
        something against a window. Same shape as the rest -- Wilson lower bound
        of p_hat over n, rescaled by the room the null leaves.

          p_hat = mass in [lo, hi] / mass in [dicermin, dicermax]
          n     = distinct_seqs
          p0    = (hi - lo + 1) / (dicermax - dicermin + 1)
        """
        size = (spec or {}).get("size")
        if not size or "lo" not in size or "hi" not in size:
            return 0.0
        counts = getattr(record, "size_counts", None) or {}
        total = sum(counts.values())
        n = getattr(record, "distinct_seqs", 0)
        if total <= 0 or n <= 0:
            return 0.0
        ext = (self.cfg.get("extractor", {}) or {})
        dmin, dmax = ext.get("dicermin", 18), ext.get("dicermax", 35)
        width = dmax - dmin + 1
        if width <= 0:
            return 0.0
        inside = sum(c for L, c in counts.items() if size["lo"] <= L <= size["hi"])
        p0 = (size["hi"] - size["lo"] + 1) / width
        return effect_strength(round((inside / total) * n), n, p0)

    def phase_cycle_strength(self, record):
        """Phasing read as the length of the run, in [0, 1].

        Same quantity as `e_phase`: cycles in register out of the window's, as a
        Wilson lower bound against the chance a position lands in register. Kept
        here rather than read off the record so it also works when the scorer is
        replayed over a .clusters.tsv, where only the columns exist.
        """
        if record.dicer_call == "N":
            return 0.0
        pl = int(record.dicer_call)
        if pl <= 0:
            return 0.0
        number = getattr(record, "phase_number", 0) or 0
        return effect_strength(min(int(number), PHASE_WINDOW_CYCLES),
                               PHASE_WINDOW_CYCLES, 1 / pl)

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
        # Size gate: every evidence term is multiplied by it, so a signal only
        # counts where the size already fits the class. Ping-pong and phasing
        # were gated from the start -- dual-strand 21-nt loci hit the 10-nt
        # overlap by geometry -- but 5'U, 5'A and the dual-strand term were not,
        # and each of them was carrying loci into a class whose window they sit
        # outside: 14.5% of Drosophila piRNA (675 of them at 18-20 nt, held by a
        # 5'U alone), 27.3% of rice siRNA, 12.1% of rice hcsiRNA. The one class
        # whose terms were all gated, phasiRNA24, sat at 0.5%.
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
            t = spec["w_1u"] * record.pct_1u * gate
            s += t; terms["1u"] = round(t, 3)
        if "w_5a" in spec:
            # 5'A -> AGO4 in plants (24-nt heterochromatic siRNA); the counterpart
            # of the 5'U/AGO1 route scored by w_1u. Unused by the arthropod config.
            t = spec["w_5a"] * getattr(record, "pct_5a", 0.0) * gate
            s += t; terms["5a"] = round(t, 3)
        if "w_fractop_mid" in spec:
            t = spec["w_fractop_mid"] * _mid(record.frac_top, sigma=sigma) * gate
            s += t; terms["fractop_mid"] = round(t, 3)
        if "w_pingpong_penalty" in spec:
            t = spec["w_pingpong_penalty"] * self._pp_signal(record) * gate
            s -= t; terms["pingpong_penalty"] = round(-t, 3)
        if "w_phase" in spec:
            # Size-gated like ping-pong: phasing only counts where the size fits,
            # which is what lets one config serve plants and arthropods.
            if spec.get("phase_on", "score") == "cycles":
                # The run in register, with no abundance in it. phase_score
                # carries a ln(mass) factor, and measured on a 24-nt class in
                # Arabidopsis -- which has no DCL5 -- it selected loci 5x deeper
                # with phase_pvalue 1.000: depth, not phasing.
                signal = self.phase_cycle_strength(record)
            else:
                # phase_score is unbounded (ratio*number*log(abun)); saturate to [0,1].
                signal = _sigmoid(record.phase_score,
                                  x0=spec.get("phase_x0", 5.0),
                                  scale=spec.get("phase_scale", 3.0))
            t = spec["w_phase"] * signal * gate
            s += t; terms["phase"] = round(t, 3)
        return s, terms

    def defining_scores(self, spec, record):
        """Evidence scores the class is defined by, in config order.

        A nested list is an OR slot: score it as the max over its members
        (e.g. piRNA is size + (5'U OR ping-pong)). A plain string is AND, as before.

        Reads the `e_*` effect strengths (in [0,1]), not the `s_*` p-value
        scores. A record without the strength attribute (e.g. the raw-feature
        stand-ins built by extras/simulate_rules.py) defaults to 0.0.
        """
        scores = []
        for name in spec.get("defining_evidence", []):
            if isinstance(name, list):
                scores.append(max(self._strength(n, spec, record) for n in name))
            else:
                scores.append(self._strength(name, spec, record))
        return scores

    def _strength(self, name, spec, record):
        """`size` is computed against the class; every other name lives on the record."""
        if name == "size":
            return self.size_strength(spec, record)
        return getattr(record, f"e_{name}", 0.0)

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
        defining = self.defining_scores(self.cfg["classes"][top_name], record)
        record.score = round(min(defining), 3) if defining else 0.0
        # No size call means no class evidence is meaningful.
        if record.dicer_call == "N":
            record.score = 0.0
        record.confidence = round(conf, 3)
        record.evidence = "; ".join(f"{k}={v}" for k, v in evidence[top_name].items())
        record.review = (
            conf < dec["review_confidence"]
            or margin < dec["min_margin"]
            or record.reads < dec["min_reads"]
            or record.dicer_call == "N"
        )
        return record
