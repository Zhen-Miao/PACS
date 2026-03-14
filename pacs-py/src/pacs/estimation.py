"""Parameter estimation with IRLS and Newton-Raphson with Firth regularization.

This module implements the core optimization algorithms for PACS, ported from
the R implementation in ``param-estimate_logit_get_p_by_t_June.R``.
"""

from __future__ import annotations

import warnings

import numpy as np
from joblib import Parallel, delayed
from scipy.special import expit

from ._types import ConvergenceStatus

# ---------------------------------------------------------------------------
# Loss, gradient, and information matrix functions
# ---------------------------------------------------------------------------

_MACHINE_EPS_SQRT = np.sqrt(np.finfo(np.float64).eps)


def loss_fun(
    p_bg: np.ndarray, q_vec: np.ndarray, y_vec: np.ndarray
) -> float:
    """Log-likelihood without Firth prior.

    Parameters
    ----------
    p_bg : np.ndarray, shape (n,)
        Open probability for each cell.
    q_vec : np.ndarray, shape (n,)
        Capturing probability for each cell.
    y_vec : np.ndarray, shape (n,)
        Observed accessibility (0 or 1) for each cell.

    Returns
    -------
    float
        Log-likelihood value.
    """
    log_pq = np.log(q_vec) + np.log(p_bg)
    log_1_pq = np.log(1.0 - q_vec * p_bg)
    return float(np.sum(y_vec * log_pq + (1.0 - y_vec) * log_1_pq))


def loss_gradient(
    X: np.ndarray,
    p_bg: np.ndarray,
    q_vec: np.ndarray,
    y_vec: np.ndarray,
) -> np.ndarray:
    """Gradient of the log-likelihood without Firth prior.

    Parameters
    ----------
    X : np.ndarray, shape (n, p)
        Design matrix.
    p_bg : np.ndarray, shape (n,)
        Open probability.
    q_vec : np.ndarray, shape (n,)
        Capturing probability.
    y_vec : np.ndarray, shape (n,)
        Observed accessibility.

    Returns
    -------
    np.ndarray, shape (p,)
        Gradient vector.
    """
    # First part: X^T (y - p)
    y_m_p = y_vec - p_bg
    grad1 = X.T @ y_m_p

    # Second part: X^T ((1 - y) * p * (1 - q) / (1 - p*q))
    fc = (1.0 - q_vec) / (1.0 - p_bg * q_vec)
    y_p = (1.0 - y_vec) * p_bg * fc
    grad2 = X.T @ y_p

    return grad1 + grad2


def infor_mat(
    X: np.ndarray, p_bg: np.ndarray, q_vec: np.ndarray
) -> np.ndarray:
    """Fisher information matrix.

    Parameters
    ----------
    X : np.ndarray, shape (n, p)
        Design matrix.
    p_bg : np.ndarray, shape (n,)
        Open probability.
    q_vec : np.ndarray, shape (n,)
        Capturing probability.

    Returns
    -------
    np.ndarray, shape (p, p)
        Information matrix.
    """
    wii = (p_bg * q_vec * (1.0 - p_bg) ** 2) / (1.0 - p_bg * q_vec)
    wii_sqrt = np.sqrt(wii)
    wii_sqrt_X = wii_sqrt[:, np.newaxis] * X
    return wii_sqrt_X.T @ wii_sqrt_X


def _compute_wii_sqrt_X(
    X: np.ndarray, p_bg: np.ndarray, q_vec: np.ndarray
) -> np.ndarray:
    """Compute sqrt(Wii) * X."""
    wii = (p_bg * q_vec * (1.0 - p_bg) ** 2) / (1.0 - p_bg * q_vec)
    wii_sqrt = np.sqrt(wii)
    return wii_sqrt[:, np.newaxis] * X


def loss_grad_pen(
    X: np.ndarray,
    p_bg: np.ndarray,
    q_vec: np.ndarray,
    inf_mat: np.ndarray,
    wii_sqrt_X: np.ndarray,
) -> np.ndarray | None:
    """Gradient of the Firth-penalized loss (invariant prior).

    Parameters
    ----------
    X : np.ndarray, shape (n, p)
        Design matrix.
    p_bg : np.ndarray, shape (n,)
        Open probability.
    q_vec : np.ndarray, shape (n,)
        Capturing probability.
    inf_mat : np.ndarray, shape (p, p)
        Information matrix.
    wii_sqrt_X : np.ndarray, shape (n, p)
        sqrt(wii) * X precomputed.

    Returns
    -------
    np.ndarray, shape (p,) or None
        Firth penalty gradient, or None if information matrix is singular.
    """
    kii = (2.0 * p_bg**2 * q_vec - 3.0 * p_bg + 1.0) / (
        1.0 - p_bg * q_vec
    )

    try:
        i_infor_m = np.linalg.solve(inf_mat, np.eye(inf_mat.shape[0]))
    except np.linalg.LinAlgError:
        return None

    # Hat matrix diagonal: h_ii = rowSums(wii_sqrt_X @ inv(I) * wii_sqrt_X)
    h_mat1 = wii_sqrt_X @ i_infor_m
    h_mat_diag = np.sum(h_mat1 * wii_sqrt_X, axis=1)

    hk = h_mat_diag * kii
    return 0.5 * (X.T @ hk)


def compute_infor_mat_tilda(
    X: np.ndarray,
    p_bg: np.ndarray,
    q_vec: np.ndarray,
    y_vec: np.ndarray,
) -> np.ndarray:
    """Alternative information matrix for Newton-Raphson.

    Parameters
    ----------
    X : np.ndarray, shape (n, p)
        Design matrix.
    p_bg : np.ndarray, shape (n,)
        Open probability.
    q_vec : np.ndarray, shape (n,)
        Capturing probability.
    y_vec : np.ndarray, shape (n,)
        Observed accessibility.

    Returns
    -------
    np.ndarray, shape (p, p)
        Alternative information matrix.
    """
    wii = (
        (-1.0 * p_bg**2 * q_vec**2)
        + q_vec * (2.0 * p_bg + y_vec - 1.0)
        - y_vec
    ) * p_bg * (1.0 - p_bg) / ((1.0 - p_bg * q_vec) ** 2)
    wii_neg = -wii
    wii_X = wii_neg[:, np.newaxis] * X
    return X.T @ wii_X


def loss_fun_star(
    X: np.ndarray,
    p_bg: np.ndarray,
    q_vec: np.ndarray,
    y_vec: np.ndarray,
) -> float:
    """Regularized log-likelihood (with Firth penalty).

    Parameters
    ----------
    X : np.ndarray, shape (n, p)
        Design matrix.
    p_bg : np.ndarray, shape (n,)
        Open probability.
    q_vec : np.ndarray, shape (n,)
        Capturing probability.
    y_vec : np.ndarray, shape (n,)
        Observed accessibility.

    Returns
    -------
    float
        Regularized log-likelihood: L + 0.5 * log(det(I)).
    """
    log_loss = loss_fun(p_bg, q_vec, y_vec)

    wii = (p_bg * q_vec * (1.0 - p_bg) ** 2) / (1.0 - p_bg * q_vec)
    wii_sqrt = np.sqrt(wii)
    wii_sqrt_X = wii_sqrt[:, np.newaxis] * X
    inf_m = wii_sqrt_X.T @ wii_sqrt_X

    sign, logdet = np.linalg.slogdet(inf_m)
    if sign <= 0:
        return float("-inf")

    return log_loss + 0.5 * logdet


# ---------------------------------------------------------------------------
# IRLS iterations
# ---------------------------------------------------------------------------


def irls_iter(
    y_vec: np.ndarray,
    X: np.ndarray,
    theta: np.ndarray,
    q_vec: np.ndarray,
    stop_criteria: float = 1e-6,
    tolerance: float = _MACHINE_EPS_SQRT,
    max_iter: int = 15,
) -> np.ndarray:
    """IRLS iteration with Firth regularization for the full model.

    Parameters
    ----------
    y_vec : np.ndarray, shape (n,)
        Observed accessibility.
    X : np.ndarray, shape (n, p)
        Design matrix.
    theta : np.ndarray, shape (p,)
        Initial parameter estimates.
    q_vec : np.ndarray, shape (n,)
        Capturing probability.
    stop_criteria : float
        Convergence threshold on sum of squared parameter changes.
    tolerance : float
        Tolerance for matrix singularity check.
    max_iter : int
        Maximum number of iterations.

    Returns
    -------
    np.ndarray, shape (p + 1,)
        Estimated parameters followed by convergence status code.
    """
    theta_est = theta.copy()
    conv_stat = ConvergenceStatus.CONVERGED
    indi = 1.0

    for _ in range(max_iter):
        p_bg = expit(X @ theta_est)

        wii_sqrt_X = _compute_wii_sqrt_X(X, p_bg, q_vec)
        inf_m = wii_sqrt_X.T @ wii_sqrt_X

        score_pen = loss_grad_pen(X, p_bg, q_vec, inf_m, wii_sqrt_X)
        if score_pen is None:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        score_total = score_pen + loss_gradient(X, p_bg, q_vec, y_vec)

        sign, logdet = np.linalg.slogdet(inf_m)
        if sign <= 0 or logdet < tolerance:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        if np.any(np.isnan(score_total)):
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        try:
            update = np.linalg.solve(inf_m, score_total)
        except np.linalg.LinAlgError:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        theta_update = theta_est + update
        indi = float(np.sum((theta_update - theta_est) ** 2))
        theta_est = theta_update

        if indi < stop_criteria:
            break
    else:
        if indi >= stop_criteria:
            conv_stat = ConvergenceStatus.MAX_ITER_REACHED

    return np.append(theta_est, float(conv_stat))


def irls_iter_null(
    y_vec: np.ndarray,
    X: np.ndarray,
    theta: np.ndarray,
    hold_zero: np.ndarray,
    q_vec: np.ndarray,
    stop_criteria: float = 1e-6,
    tolerance: float = _MACHINE_EPS_SQRT,
    max_iter: int = 15,
) -> np.ndarray:
    """IRLS iteration with Firth regularization for the null model.

    Parameters are held at zero for indices in ``hold_zero``.

    Parameters
    ----------
    y_vec : np.ndarray, shape (n,)
        Observed accessibility.
    X : np.ndarray, shape (n, p)
        Design matrix (full, not reduced).
    theta : np.ndarray, shape (p,)
        Initial parameter estimates.
    hold_zero : np.ndarray
        Indices of parameters to hold at zero (0-based).
    q_vec : np.ndarray, shape (n,)
        Capturing probability.
    stop_criteria : float
        Convergence threshold.
    tolerance : float
        Tolerance for singularity check.
    max_iter : int
        Maximum iterations.

    Returns
    -------
    np.ndarray, shape (p + 1,)
        Estimated parameters followed by convergence status code.
    """
    theta_est = theta.copy()
    conv_stat = ConvergenceStatus.CONVERGED
    indi = 1.0

    # Indices of parameters that are free to vary
    all_idx = np.arange(len(theta))
    poi = np.setdiff1d(all_idx, hold_zero)

    for _ in range(max_iter):
        p_bg = expit(X @ theta_est)

        wii_sqrt_X = _compute_wii_sqrt_X(X, p_bg, q_vec)
        inf_m = wii_sqrt_X.T @ wii_sqrt_X

        sign, logdet = np.linalg.slogdet(inf_m)
        if sign <= 0 or logdet < tolerance:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        score_pen = loss_grad_pen(X, p_bg, q_vec, inf_m, wii_sqrt_X)
        if score_pen is None:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        score_total = score_pen + loss_gradient(X, p_bg, q_vec, y_vec)

        if np.any(np.isnan(score_total)):
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        # Only invert the submatrix for free parameters
        inf_sub = inf_m[np.ix_(poi, poi)]
        try:
            update_free = np.linalg.solve(inf_sub, score_total[poi])
        except np.linalg.LinAlgError:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        full_update = np.zeros_like(theta_est)
        full_update[poi] = update_free

        theta_update = theta_est + full_update
        indi = float(np.sum((theta_update - theta_est) ** 2))
        theta_est = theta_update

        if indi < stop_criteria:
            break
    else:
        if indi >= stop_criteria:
            conv_stat = ConvergenceStatus.MAX_ITER_REACHED

    return np.append(theta_est, float(conv_stat))


def irls_iter_nt(
    y_vec: np.ndarray,
    X: np.ndarray,
    theta: np.ndarray,
    q_vec: np.ndarray,
    stop_criteria: float = 1e-5,
    tolerance: float = _MACHINE_EPS_SQRT,
    max_iter: int = 15,
) -> np.ndarray:
    """Newton-Raphson iteration with Firth regularization (full model).

    Uses a blended information matrix (0.5 * I_tilda + 0.5 * I) and monitors
    loss for convergence.

    Parameters
    ----------
    y_vec : np.ndarray, shape (n,)
        Observed accessibility.
    X : np.ndarray, shape (n, p)
        Design matrix.
    theta : np.ndarray, shape (p,)
        Initial parameter estimates.
    q_vec : np.ndarray, shape (n,)
        Capturing probability.
    stop_criteria : float
        Convergence threshold.
    tolerance : float
        Tolerance for singularity check.
    max_iter : int
        Maximum iterations.

    Returns
    -------
    np.ndarray, shape (p + 1,)
        Estimated parameters followed by convergence status code.
    """
    theta_est = theta.copy()
    conv_stat = ConvergenceStatus.CONVERGED
    loss_star_hist = np.full(max_iter, np.nan)

    for n_iter in range(max_iter):
        p_bg = expit(X @ theta_est)

        inf_mat_tilda = compute_infor_mat_tilda(X, p_bg, q_vec, y_vec)
        wii_sqrt_X = _compute_wii_sqrt_X(X, p_bg, q_vec)
        inf_m = wii_sqrt_X.T @ wii_sqrt_X

        loss_star_hist[n_iter] = loss_fun_star(X, p_bg, q_vec, y_vec)

        # Blend information matrices
        inf_mat_tilda = 0.5 * inf_mat_tilda + 0.5 * inf_m

        sign, logdet = np.linalg.slogdet(inf_mat_tilda)
        if sign <= 0 or logdet < tolerance:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        score_pen = loss_grad_pen(X, p_bg, q_vec, inf_mat_tilda, wii_sqrt_X)
        if score_pen is None:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        score_total = score_pen + loss_gradient(X, p_bg, q_vec, y_vec)

        if np.any(np.isnan(score_total)):
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        try:
            update = np.linalg.solve(inf_mat_tilda, score_total)
        except np.linalg.LinAlgError:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        theta_update = theta_est + update

        # Monitor loss for convergence
        if n_iter >= 1 and not np.isnan(loss_star_hist[n_iter]) and loss_star_hist[n_iter] != float("-inf"):
            diff = loss_star_hist[n_iter] - loss_star_hist[n_iter - 1]
            if diff > 0:
                theta_est = theta_update
                if diff < 1.0:
                    break
            else:
                conv_stat = ConvergenceStatus.LOSS_NOT_INCREASING
                break
        elif not np.isnan(loss_star_hist[n_iter]) and loss_star_hist[n_iter] != float("-inf"):
            theta_est = theta_update
        else:
            conv_stat = ConvergenceStatus.LOSS_NOT_INCREASING
            break
    else:
        conv_stat = ConvergenceStatus.MAX_ITER_REACHED

    return np.append(theta_est, float(conv_stat))


def irls_iter_nt_null(
    y_vec: np.ndarray,
    X: np.ndarray,
    theta: np.ndarray,
    hold_zero: np.ndarray,
    q_vec: np.ndarray,
    stop_criteria: float = 1e-5,
    tolerance: float = _MACHINE_EPS_SQRT,
    max_iter: int = 15,
) -> np.ndarray:
    """Newton-Raphson iteration with Firth regularization (null model).

    Parameters at ``hold_zero`` indices are constrained to zero.

    Parameters
    ----------
    y_vec : np.ndarray, shape (n,)
        Observed accessibility.
    X : np.ndarray, shape (n, p)
        Design matrix.
    theta : np.ndarray, shape (p,)
        Initial parameter estimates.
    hold_zero : np.ndarray
        Indices of parameters to hold at zero (0-based).
    q_vec : np.ndarray, shape (n,)
        Capturing probability.
    stop_criteria : float
        Convergence threshold.
    tolerance : float
        Tolerance for singularity check.
    max_iter : int
        Maximum iterations.

    Returns
    -------
    np.ndarray, shape (p + 1,)
        Estimated parameters followed by convergence status code.
    """
    theta_est = theta.copy()
    conv_stat = ConvergenceStatus.CONVERGED
    loss_star_hist = np.full(max_iter, np.nan)

    all_idx = np.arange(len(theta))
    poi = np.setdiff1d(all_idx, hold_zero)

    for n_iter in range(max_iter):
        p_bg = expit(X @ theta_est)

        inf_mat_tilda = compute_infor_mat_tilda(X, p_bg, q_vec, y_vec)
        wii_sqrt_X = _compute_wii_sqrt_X(X, p_bg, q_vec)
        inf_m = wii_sqrt_X.T @ wii_sqrt_X

        loss_star_hist[n_iter] = loss_fun_star(X, p_bg, q_vec, y_vec)

        # Blend
        inf_mat_tilda = 0.5 * inf_mat_tilda + 0.5 * inf_m

        sign, logdet = np.linalg.slogdet(inf_mat_tilda)
        if sign <= 0 or logdet < tolerance:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        score_pen = loss_grad_pen(X, p_bg, q_vec, inf_mat_tilda, wii_sqrt_X)
        if score_pen is None:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        score_total = score_pen + loss_gradient(X, p_bg, q_vec, y_vec)

        if np.any(np.isnan(score_total)):
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        inf_sub = inf_mat_tilda[np.ix_(poi, poi)]
        try:
            update_free = np.linalg.solve(inf_sub, score_total[poi])
        except np.linalg.LinAlgError:
            conv_stat = ConvergenceStatus.SINGULAR_MATRIX
            break

        full_update = np.zeros_like(theta_est)
        full_update[poi] = update_free
        theta_update = theta_est + full_update

        # Monitor loss
        if n_iter >= 1 and not np.isnan(loss_star_hist[n_iter]) and loss_star_hist[n_iter] != float("-inf"):
            diff = loss_star_hist[n_iter] - loss_star_hist[n_iter - 1]
            if diff > 0:
                theta_est = theta_update
                if diff < 1.0:
                    break
            else:
                conv_stat = ConvergenceStatus.LOSS_NOT_INCREASING
                break
        elif not np.isnan(loss_star_hist[n_iter]) and loss_star_hist[n_iter] != float("-inf"):
            theta_est = theta_update
        else:
            conv_stat = ConvergenceStatus.LOSS_NOT_INCREASING
            break
    else:
        conv_stat = ConvergenceStatus.MAX_ITER_REACHED

    return np.append(theta_est, float(conv_stat))


# ---------------------------------------------------------------------------
# Orchestration: estimate_parameters / estimate_parameters_null
# ---------------------------------------------------------------------------


def _run_irls_single(y_vec, X, par_initial, q_vec):
    """Run IRLS for a single feature."""
    return irls_iter(y_vec, X, par_initial, q_vec)


def _run_irls_null_single(y_vec, X, par_initial, hold_zero, q_vec):
    """Run IRLS null for a single feature."""
    return irls_iter_null(y_vec, X, par_initial, hold_zero, q_vec)


def _run_nt_single(y_vec, X, par_initial, q_vec):
    """Run Newton-Raphson for a single feature."""
    return irls_iter_nt(y_vec, X, par_initial, q_vec)


def _run_nt_null_single(y_vec, X, par_initial, hold_zero, q_vec):
    """Run Newton-Raphson null for a single feature."""
    return irls_iter_nt_null(y_vec, X, par_initial, hold_zero, q_vec)


def estimate_parameters(
    r_by_c: np.ndarray,
    design_mat: np.ndarray,
    par_initial: np.ndarray,
    cap_rate_vec: np.ndarray,
    n_jobs: int = 1,
    verbose: bool = True,
) -> np.ndarray:
    """Estimate parameters for the full model.

    Runs IRLS with Firth regularization for each feature in parallel.
    Features that fail to converge are retried with Newton-Raphson.

    Parameters
    ----------
    r_by_c : np.ndarray, shape (n_features, n_cells)
        Region-by-cell count matrix.
    design_mat : np.ndarray, shape (n_cells, n_params)
        Design matrix.
    par_initial : np.ndarray, shape (n_params,)
        Initial parameter values.
    cap_rate_vec : np.ndarray, shape (n_cells,)
        Capturing probability for each cell.
    n_jobs : int
        Number of parallel jobs.
    verbose : bool
        Whether to print progress.

    Returns
    -------
    np.ndarray, shape (n_params + 1, n_features)
        Rows are parameters followed by convergence status.
        Columns are features.
    """
    n_features = r_by_c.shape[0]
    # Transpose to cell-by-region
    c_by_r = r_by_c.T

    # Run IRLS for each feature
    results = Parallel(n_jobs=n_jobs)(
        delayed(_run_irls_single)(
            c_by_r[:, i], design_mat, par_initial, cap_rate_vec
        )
        for i in range(n_features)
    )

    theta_summaries = np.column_stack(results)

    # Retry failed features with Newton-Raphson
    conv_row = theta_summaries[-1, :]
    fail_mask = conv_row != ConvergenceStatus.CONVERGED
    n_fail = int(np.sum(fail_mask))

    if verbose and n_fail > 0:
        print(f"{n_fail} regions failed to converge, trying Newton-Raphson")

    if n_fail > 0:
        fail_indices = np.where(fail_mask)[0]
        par_neg = -par_initial

        nt_results = Parallel(n_jobs=1)(
            delayed(_run_nt_single)(
                c_by_r[:, i], design_mat, par_neg, cap_rate_vec
            )
            for i in fail_indices
        )

        for j, idx in enumerate(fail_indices):
            theta_summaries[:, idx] = nt_results[j]

        if verbose:
            still_fail = sum(
                1 for r in nt_results
                if r[-1] != ConvergenceStatus.CONVERGED
            )
            if still_fail > 0:
                print(
                    f"{still_fail} regions failed to converge after Newton-Raphson"
                )

    return theta_summaries


def estimate_parameters_null(
    r_by_c: np.ndarray,
    design_mat: np.ndarray,
    par_initial: np.ndarray,
    hold_zero: np.ndarray,
    cap_rate_vec: np.ndarray,
    n_jobs: int = 1,
    verbose: bool = True,
) -> np.ndarray:
    """Estimate parameters for the null model.

    Parameters at ``hold_zero`` indices are constrained to zero.

    Parameters
    ----------
    r_by_c : np.ndarray, shape (n_features, n_cells)
        Region-by-cell count matrix.
    design_mat : np.ndarray, shape (n_cells, n_params)
        Design matrix (full, not reduced).
    par_initial : np.ndarray, shape (n_params,)
        Initial parameter values (parameters at hold_zero should be 0).
    hold_zero : np.ndarray
        Indices of parameters to hold at zero (0-based).
    cap_rate_vec : np.ndarray, shape (n_cells,)
        Capturing probability for each cell.
    n_jobs : int
        Number of parallel jobs.
    verbose : bool
        Whether to print progress.

    Returns
    -------
    np.ndarray, shape (n_params + 1, n_features)
        Rows are parameters followed by convergence status.
    """
    n_features = r_by_c.shape[0]
    c_by_r = r_by_c.T

    results = Parallel(n_jobs=n_jobs)(
        delayed(_run_irls_null_single)(
            c_by_r[:, i], design_mat, par_initial, hold_zero, cap_rate_vec
        )
        for i in range(n_features)
    )

    theta_summaries = np.column_stack(results)

    # Retry failures with Newton-Raphson
    conv_row = theta_summaries[-1, :]
    fail_mask = conv_row != ConvergenceStatus.CONVERGED
    n_fail = int(np.sum(fail_mask))

    if verbose and n_fail > 0:
        print(f"{n_fail} regions failed to converge, trying Newton-Raphson")

    if n_fail > 0:
        fail_indices = np.where(fail_mask)[0]
        par_neg = -par_initial

        nt_results = Parallel(n_jobs=1)(
            delayed(_run_nt_null_single)(
                c_by_r[:, i], design_mat, par_neg, hold_zero, cap_rate_vec
            )
            for i in fail_indices
        )

        for j, idx in enumerate(fail_indices):
            theta_summaries[:, idx] = nt_results[j]

        if verbose:
            still_fail = sum(
                1 for r in nt_results
                if r[-1] != ConvergenceStatus.CONVERGED
            )
            if still_fail > 0:
                print(
                    f"{still_fail} regions failed to converge after Newton-Raphson"
                )

    return theta_summaries
