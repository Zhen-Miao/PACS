"""Differential identification via generalized Likelihood Ratio Test.

Ported from ``differential_identification.R``.
"""

from __future__ import annotations

import numpy as np
from joblib import Parallel, delayed
from scipy.special import expit
from scipy.stats import chi2


def loss_firth_from_wii(wii_sqrt: np.ndarray, X: np.ndarray) -> float:
    """Compute the Firth regularization term from sqrt(Wii) weights.

    Parameters
    ----------
    wii_sqrt : np.ndarray, shape (n,)
        Square root of the diagonal weight elements.
    X : np.ndarray, shape (n, p)
        Design matrix.

    Returns
    -------
    float
        0.5 * log(det(X^T W X)).
    """
    wii_sqrt_X = wii_sqrt[:, np.newaxis] * X
    inf_mat = wii_sqrt_X.T @ wii_sqrt_X
    sign, logdet = np.linalg.slogdet(inf_mat)
    if sign <= 0:
        return float("-inf")
    return 0.5 * logdet


def compare_models(
    X_full: np.ndarray,
    theta_full: np.ndarray,
    X_null: np.ndarray,
    theta_null: np.ndarray,
    q_vec: np.ndarray,
    c_by_r: np.ndarray,
    df_test: int,
    n_jobs: int = 1,
) -> np.ndarray:
    """Use gLRT to compare full and null models and obtain p-values.

    Parameters
    ----------
    X_full : np.ndarray, shape (n_cells, p_full)
        Design matrix for the full model.
    theta_full : np.ndarray, shape (p_full, n_features)
        Estimated coefficients for the full model.
    X_null : np.ndarray, shape (n_cells, p_null)
        Design matrix for the null model.
    theta_null : np.ndarray, shape (p_null, n_features)
        Estimated coefficients for the null model.
    q_vec : np.ndarray, shape (n_cells,)
        Capturing probability for each cell.
    c_by_r : np.ndarray, shape (n_cells, n_features)
        Cell-by-region count matrix.
    df_test : int
        Degrees of freedom for the test.
    n_jobs : int
        Number of parallel jobs.

    Returns
    -------
    np.ndarray, shape (n_features,)
        P-values from the generalized likelihood ratio test.
    """
    n_features = theta_full.shape[1]

    # --- Full model ---
    # Vectorized across all features: X @ theta -> (n_cells, n_features)
    x_times_theta = X_full @ theta_full
    p_bg = expit(x_times_theta)

    # Log-likelihood without Firth (vectorized)
    log_pq = np.log(q_vec)[:, np.newaxis] + np.log(p_bg)
    log_1_pq = np.log(1.0 - p_bg * q_vec[:, np.newaxis])
    log_loss_full = np.sum(
        c_by_r * log_pq + (1.0 - c_by_r) * log_1_pq, axis=0
    )

    # Weights for Firth penalty
    wii = (
        p_bg * (1.0 - p_bg) ** 2 * q_vec[:, np.newaxis]
    ) / (1.0 - p_bg * q_vec[:, np.newaxis])
    wii_sqrt_full = np.sqrt(wii)

    # --- Null model ---
    x_times_theta_null = X_null @ theta_null
    p_bg_null = expit(x_times_theta_null)

    log_pq_null = np.log(q_vec)[:, np.newaxis] + np.log(p_bg_null)
    log_1_pq_null = np.log(1.0 - p_bg_null * q_vec[:, np.newaxis])
    log_loss_null = np.sum(
        c_by_r * log_pq_null + (1.0 - c_by_r) * log_1_pq_null, axis=0
    )

    wii_null = (
        p_bg_null * (1.0 - p_bg_null) ** 2 * q_vec[:, np.newaxis]
    ) / (1.0 - p_bg_null * q_vec[:, np.newaxis])
    wii_sqrt_null = np.sqrt(wii_null)

    # Compute Firth penalty per feature (requires per-feature matrix ops)
    firth_full = Parallel(n_jobs=n_jobs)(
        delayed(loss_firth_from_wii)(wii_sqrt_full[:, i], X_full)
        for i in range(n_features)
    )
    firth_null = Parallel(n_jobs=n_jobs)(
        delayed(loss_firth_from_wii)(wii_sqrt_null[:, i], X_null)
        for i in range(n_features)
    )

    log_loss_star_full = np.array(firth_full) + log_loss_full
    log_loss_star_null = np.array(firth_null) + log_loss_null

    # Test statistic and p-value
    test_stat = 2.0 * (log_loss_star_full - log_loss_star_null)
    p_values = chi2.sf(test_stat, df=df_test)

    return p_values
