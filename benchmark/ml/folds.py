"""Chromosome-grouped leave-one-unit-out folds."""
from typing import Iterator

import numpy as np
import pandas as pd

from benchmark.ml.dataset import load_ml_config


def chromosome_units(table: pd.DataFrame, min_share: float = 0.05) -> pd.Series:
    """Fold unit per row: "<organism>:<chrom>" when the chromosome holds >= min_share of its
    organism's loci, else "<organism>:other". Shares are computed per organism."""
    org = table["organism"].astype(str)
    chrom = table["chrom"].astype(str)
    share = chrom.groupby(org).transform(lambda s: s.map(s.value_counts(normalize=True)))
    big = share >= min_share
    return (org + ":" + chrom.where(big, "other")).rename("unit")


def split(table: pd.DataFrame, min_share: float | None = None
          ) -> Iterator[tuple[np.ndarray, np.ndarray, str]]:
    """Leave-one-unit-out within each profile; yields (train_idx, val_idx, unit) as index
    labels of `table`. min_share defaults to folds.min_share of benchmark/ml/config.yaml.
    The table is assumed to hold tuning libraries only (build_table refuses test ones)."""
    if min_share is None:
        min_share = float(load_ml_config()["folds"]["min_share"])
    for _, sub in table.groupby("profile", sort=True):
        units = chromosome_units(sub, min_share)
        names = sorted(units.unique())
        if len(names) < 2:
            continue
        for unit in names:
            val = units == unit
            yield sub.index[~val.to_numpy()].to_numpy(), sub.index[val.to_numpy()].to_numpy(), unit
