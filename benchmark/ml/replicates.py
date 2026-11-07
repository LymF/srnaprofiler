"""Replicate reproducibility features: how often, and how strongly, a locus recurs."""
import numpy as np
import pandas as pd

from benchmark.metrics_class import IntervalIndex
from benchmark.ml.controls import overlaps

COLS = ["repro_frac", "repro_rpm_ratio"]


def repro_features(table: pd.DataFrame, min_overlap_frac: float | None = None) -> pd.DataFrame:
    """repro_frac: fraction of the other libraries of the (organism, tissue) group holding a locus
    that overlaps this one (match_loci rule). repro_rpm_ratio: median RPM of the matched loci over
    this locus's RPM (NaN without a match or with RPM <= 0). NaN for both with a single library."""
    if min_overlap_frac is None:
        from benchmark.ml.dataset import load_ml_config
        min_overlap_frac = float(load_ml_config()["min_overlap_frac"])
    out = pd.DataFrame(np.nan, index=table.index, columns=COLS)
    for _, grp in table.groupby(["organism", "tissue"], sort=False):
        accs = list(dict.fromkeys(grp["acc"]))
        if len(accs) < 2:
            continue
        per_acc = {}
        for a in accs:
            g = grp[grp["acc"] == a]
            iv = [(c, int(s), int(e)) for c, s, e in zip(g["chrom"], g["start"], g["end"])]
            per_acc[a] = (IntervalIndex(iv) if iv else None, iv, g["reads_rpm"].to_numpy(float))
        for a in accs:
            for row in grp.index[grp["acc"] == a]:
                c, s, e = table.at[row, "chrom"], int(table.at[row, "start"]), int(table.at[row, "end"])
                n_hit, rpms = 0, []
                for b in accs:
                    if b == a or per_acc[b][0] is None:
                        continue
                    index, iv, rpm = per_acc[b]
                    hits = list(overlaps(index, iv, c, s, e, min_overlap_frac))
                    if hits:
                        n_hit += 1
                        rpms.extend(rpm[hits])
                out.at[row, "repro_frac"] = n_hit / (len(accs) - 1)
                own = table.at[row, "reads_rpm"]
                if rpms and own and own > 0:
                    out.at[row, "repro_rpm_ratio"] = float(np.nanmedian(rpms)) / own
    return out
