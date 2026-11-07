"""Statistics for the model comparison."""
import math

import numpy as np
from scipy import stats as _st


def corrected_resampled_ttest(a: list[float], b: list[float], n_train: int, n_test: int
                              ) -> tuple[float, float]:
    """Nadeau & Bengio (2003) corrected resampled t-test on paired scores.

    t = mean(d) / sqrt((1/k + n_test/n_train) * var(d)), d = a - b over k folds, var with
    ddof=1, two-sided p from Student t with k-1 degrees of freedom. Zero variance gives
    (0, 1) when the mean difference is 0 and (+-inf, 0) otherwise; k < 2 gives (nan, nan).
    """
    d = np.asarray(a, float) - np.asarray(b, float)
    k = len(d)
    if k < 2:
        return float("nan"), float("nan")
    mean, var = float(d.mean()), float(d.var(ddof=1))
    if var <= 1e-24:
        return (0.0, 1.0) if abs(mean) <= 1e-12 else (math.copysign(math.inf, mean), 0.0)
    t = mean / math.sqrt((1.0 / k + n_test / n_train) * var)
    return t, float(2 * _st.t.sf(abs(t), k - 1))
