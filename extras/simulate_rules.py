#!/usr/bin/env python3
"""Re-score an existing .clusters.tsv under a different rules.yaml.

The TSV already holds every feature the decider reads, so a rule change can be
evaluated by replaying the scorer over it -- no BAM, no re-measuring. A full run
takes minutes; this takes seconds. Replaying the run's own config reproduces its
classes and review flags exactly, which is what makes the comparison meaningful.

Usage:
    python3 extras/simulate_rules.py run.clusters.tsv candidate.yaml [baseline.yaml]

Without `baseline.yaml` the TSV's own columns are the baseline.
"""
import sys
import types

import pandas as pd

sys.path.insert(0, __file__.rsplit("/", 2)[0])
from decider import RuleScorer, load_config  # noqa: E402

FEATURES = ("dicer_call", "frac_top", "pct_1u", "pct_5a", "pingpong_frac",
            "phase_score", "reads")


def replay(df, config):
    """Run the scorer over every row; returns (classes, review flags)."""
    scorer = RuleScorer(config)
    cls, review = [], []
    missing = [f for f in FEATURES if f not in df.columns]
    if missing:
        # Silently defaulting a missing feature to 0 makes the replay disagree
        # with the run for reasons that have nothing to do with the rule change.
        sys.exit(f"{', '.join(missing)} missing from the TSV — it predates these "
                 f"columns; re-run classify before simulating.")
    for row in df.itertuples(index=False):
        r = types.SimpleNamespace(**{f: getattr(row, f) for f in FEATURES})
        r.dicer_call = str(r.dicer_call)
        scorer.decide(r)
        cls.append(r.cls)
        review.append(r.review)
    return pd.Series(cls, index=df.index), pd.Series(review, index=df.index)


def main(tsv_path, candidate_path, baseline_path=None):
    d = pd.read_csv(tsv_path, sep="\t", low_memory=False)
    if baseline_path:
        base_cls, base_review = replay(d, load_config(baseline_path))
    else:
        base_cls = d.cls
        base_review = d.review.astype(str).str.lower() == "true"
    new_cls, new_review = replay(d, load_config(candidate_path))

    changed = base_cls != new_cls
    print(f"{len(d)} loci; {changed.sum()} change class ({changed.mean() * 100:.2f}%)")
    print(f"review: {base_review.mean():.3f} -> {new_review.mean():.3f}\n")

    print("class shares:")
    for cls in sorted(set(base_cls) | set(new_cls)):
        print(f"   {cls:8s} {(base_cls == cls).mean():.3f} -> {(new_cls == cls).mean():.3f}")

    if changed.any():
        print("\nmoves (baseline -> candidate), with the read mass involved:")
        moves = pd.crosstab(base_cls[changed], new_cls[changed])
        for a in moves.index:
            for b in moves.columns:
                if moves.loc[a, b]:
                    mass = d.reads[changed & (base_cls == a) & (new_cls == b)].sum()
                    print(f"   {a} -> {b}: {moves.loc[a, b]} loci, {mass:.0f} reads")

        dc = pd.to_numeric(d.dicer_call, errors="coerce")
        print("\nsize of the loci that moved:")
        print("   " + ", ".join(f"{int(k)}nt: {v}"
                                for k, v in dc[changed].value_counts().head(8).items()))


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(*sys.argv[1:4])
