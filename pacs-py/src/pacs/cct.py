"""Cauchy Combination Test for p-value aggregation.

Adopted from Liu and Xie, JASA 2020. Please cite their publication if you
use this function.

Ported from ``CCT_function.R``.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import cauchy


def cauchy_combination_test(p_values: np.ndarray) -> np.ndarray:
    """Cauchy Combination Test for aggregating p-values across tests.

    Parameters
    ----------
    p_values : np.ndarray, shape (n_features, n_tests)
        Matrix of p-values where each row is a feature and each column
        is a p-value from an individual test.

    Returns
    -------
    np.ndarray, shape (n_features,)
        Combined p-value for each feature.
    """
    p_values = np.atleast_2d(np.asarray(p_values, dtype=np.float64))
    n_features, n_tests = p_values.shape

    weights = np.full_like(p_values, 1.0 / n_tests)

    # Identify rows with very small p-values
    has_small = np.any(p_values < 1e-16, axis=1)

    cct_stat = np.zeros(n_features)

    # Rows without very small p-values: standard formula
    if np.any(~has_small):
        normal_rows = ~has_small
        cct_stat[normal_rows] = np.sum(
            weights[normal_rows] * np.tan((0.5 - p_values[normal_rows]) * np.pi),
            axis=1,
        )

    # Rows with very small p-values: use asymptotic approximation
    if np.any(has_small):
        for i in np.where(has_small)[0]:
            small_mask = p_values[i] < 1e-16
            # For very small p-values: w/p / pi
            stat = np.sum(
                weights[i, small_mask] / p_values[i, small_mask]
            ) / np.pi
            # For normal p-values: standard tan formula
            stat += np.sum(
                weights[i, ~small_mask]
                * np.tan((0.5 - p_values[i, ~small_mask]) * np.pi)
            )
            cct_stat[i] = stat

    # Convert to p-values
    is_large = cct_stat > 1e15
    result = 1.0 - cauchy.cdf(cct_stat)

    # For very large statistics, use tail approximation
    if np.any(is_large):
        result[is_large] = 1.0 / (cct_stat[is_large] * np.pi)

    return result
