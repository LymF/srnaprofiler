"""Pluggable decider. v1: RuleScorer (transparent scoring model)."""
import math
from collections.abc import Mapping

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
    # Reported instead of a class when the size fits no class's window. Distinct
    # from `review`, which means the classes compete: here none of them applies.
    AMBIGUOUS = "ambiguous"

    @staticmethod
    def enabled_classes(config):
        """Classes the config leaves switched on, in order.

        A class is a description of a pathway, and a pathway can be absent from
        a clade: 24-nt reproductive phasiRNA needs DCL5, which monocots have and
        Arabidopsis does not. `enabled: false` removes the class everywhere it
        would otherwise act -- scoring here, and island merging through
        `classify.class_size_windows`.
        """
        return {name: spec for name, spec in (config.get("classes") or {}).items()
                if (spec or {}).get("enabled", True)}

    def __init__(self, config):
        self.cfg = config
        # Configurable shape parameters (defaults = previous behavior).
        self.shape = (config.get("decision", {}) or {}).get("shape", {}) or {}
        if not self.enabled_classes(config):
            raise ValueError("no enabled class in the config: every locus would "
                             "come out ambiguous")
        for name, spec in self.enabled_classes(config).items():
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

    @staticmethod
    def _phase_supported(spec, record):
        """Whether the record's phasing passes the class's `phase_max_pvalue`.

        Unset: always. Set: a missing, NaN or larger p-value means no phasing
        signal. `phase_number` alone counts runs in register that chance makes
        too, and in a class with a deep, abundant background (24-nt hcsiRNA)
        chance is plenty.
        """
        limit = spec.get("phase_max_pvalue")
        if limit is None:
            return True
        try:
            p = float(getattr(record, "phase_pvalue", None))
        except (TypeError, ValueError):
            return False
        return not math.isnan(p) and p <= limit

    def _score_class(self, name, spec, record):
        size = self._size(record)
        edge = self.shape.get("size_edge", 3.0)
        sigma = self.shape.get("fractop_sigma", 0.15)
        # Size gate: every evidence term is multiplied by it, so a signature only
        # counts where the size already fits the class.
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
            if not self._phase_supported(spec, record):
                signal = 0.0
            elif spec.get("phase_on", "score") == "cycles":
                # The run in register, with no abundance in it: phase_score's ln(mass)
                # factor selects depth in a class that dominates the genome.
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
        """Evidence scores the class is defined by, in config order. A nested list is
        an OR slot, scored as the max over its members; a plain string is AND.
        Reads the `e_*` effect strengths, not the `s_*` p-values.
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
        if name == "phase" and not self._phase_supported(spec, record):
            return 0.0
        return getattr(record, f"e_{name}", 0.0)

    def decide(self, record):
        scores = {}
        evidence = {}
        classes = self.enabled_classes(self.cfg)
        # A class whose gated defining evidence fails does not compete: scored,
        # it still enters the softmax through its size term and splits the
        # confidence of every locus of a neighbouring class. If that leaves
        # nothing, score them all (defensive; the gate then only zeroes terms).
        competing = {n: sp for n, sp in classes.items()
                     if self._phase_supported(sp, record)} or classes
        for name, spec in competing.items():
            scores[name], evidence[name] = self._score_class(name, spec, record)
        # Excluded classes are reported as None (NA in the .clusters.tsv).
        record.class_scores = {k: (round(scores[k], 3) if k in scores else None)
                               for k in classes}
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
        if self._outside_window(self.cfg["classes"][top_name], record):
            # No class window contains this locus. With every term size-gated an "N"
            # scores zero everywhere, so the argmax would be the config's first key.
            # Report the size verdict instead, keeping the nearest class in the trail.
            record.cls = self.AMBIGUOUS
            record.score = None
            # review stays True: a locus whose size fits no class is exactly what it
            # means, and clearing it would let these pass a `review == False` filter.
            record.review = True
            record.evidence = f"nearest={top_name}; " + record.evidence
        else:
            failed = self._cluster_failures(self.cfg["classes"][top_name], record)
            if failed:
                # Same downgrade as above: the class still names the nearest, but
                # the locus is not a cluster of it.
                record.cls = self.AMBIGUOUS
                record.score = None
                record.review = True
                record.evidence = (f"nearest={top_name}; "
                                   f"below_cluster_criterion({','.join(failed)}); "
                                   + record.evidence)
        return record

    @staticmethod
    def _cluster_failures(spec, record):
        """Which of a class's optional cluster criteria the locus fails.

        `cluster_min_length` (bp, on the span of the reads, `read_span`; falls back to
        the window end - start when the record has none) and `cluster_max_density_pvalue` (proTRAC,
        Rosenkranz & Zischler 2012). Either key absent: not applied. A missing or
        NaN density (e.g. a replayed .clusters.tsv without the column) leaves
        only the length part.
        """
        failed = []
        min_len = spec.get("cluster_min_length")
        if min_len is not None:
            span = getattr(record, "read_span", None)
            if span is None:   # replayed .clusters.tsv without the column
                span = record.end - record.start
            if span < min_len:
                failed.append("length")
        max_p = spec.get("cluster_max_density_pvalue")
        if max_p is not None:
            try:
                p = float(getattr(record, "density_pvalue", None))
            except (TypeError, ValueError):
                p = float("nan")
            if not math.isnan(p) and p > max_p:
                failed.append("density")
        return failed

    def _outside_window(self, spec, record):
        """True when the locus's dominant size is not inside this class's window."""
        window = spec.get("size")
        if not window:
            return False
        size = self._size(record)
        return size is None or not (window["lo"] <= size <= window["hi"])

    def passes_gate(self, name, record):
        """Whether a locus shows the defining biology of class `name` (spec v2 s2.1).

        The checks the rules apply to their own winner, from the same config:
        size window, cluster length/density, phasing p-value, and the class's
        defining evidence at or above `decision.ml_gate_min_score`. No size call
        ("N") passes nothing. Pure: the record is not touched.
        """
        spec = self.enabled_classes(self.cfg).get(name)
        if spec is None or record.dicer_call == "N":
            return False
        if self._outside_window(spec, record):
            return False
        if self._cluster_failures(spec, record):
            return False
        if not self._phase_supported(spec, record):
            return False
        if spec.get("defining_evidence"):
            floor = (self.cfg.get("decision", {}) or {}).get("ml_gate_min_score", 0.5)
            if min(self.defining_scores(spec, record)) < floor:
                return False
        return True

    def gated_choice(self, proba: Mapping, record):
        """Highest-probability class the locus is allowed to have, else AMBIGUOUS.

        `other` on top is the model abstaining; miRNA comes only from miRDeep2,
        so it is skipped; a class failing the gate yields to the next one.
        """
        for name, _ in sorted(proba.items(), key=lambda kv: kv[1], reverse=True):
            if name == "other":
                return self.AMBIGUOUS
            if name == "miRNA":
                continue
            if self.passes_gate(name, record):
                return name
        return self.AMBIGUOUS
