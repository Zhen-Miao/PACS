"""Main PACS testing functions for differential accessible region identification.

Ported from ``PACS_test_auto.R``, ``PACS_test_logit.R``, and
``PACS_test_cumulative.R``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._types import PACSResult
from ._utils import build_design_matrix
from .differential import compare_models
from .estimation import estimate_parameters, estimate_parameters_null


def pacs_test_logit(
    metadata: pd.DataFrame,
    formula_full: str,
    formula_null: str,
    pic_matrix: np.ndarray,
    cap_rates: np.ndarray,
    n_jobs: int = 1,
    par_initial_null: np.ndarray | None = None,
    par_initial_full: np.ndarray | None = None,
) -> PACSResult:
    """Profile likelihood ratio test using binary logit model.

    Parameters
    ----------
    metadata : pd.DataFrame
        Cell metadata with columns referenced in the formulas.
    formula_full : str
        R-style formula for the full model, e.g., ``'~ cell_type + batch'``.
    formula_null : str
        R-style formula for the null model, e.g., ``'~ batch'``.
    pic_matrix : np.ndarray, shape (n_features, n_cells)
        Binary region-by-cell matrix.
    cap_rates : np.ndarray, shape (n_cells,)
        Capturing probability for each cell.
    n_jobs : int
        Number of parallel jobs.
    par_initial_null : np.ndarray, optional
        Initial parameter values for the null model.
    par_initial_full : np.ndarray, optional
        Initial parameter values for the full model.

    Returns
    -------
    PACSResult
        Test results with p-values and convergence status.
    """
    # Build design matrices
    X_full, col_full = build_design_matrix(formula_full, metadata)
    X_null, col_null = build_design_matrix(formula_null, metadata)

    # Reorder full matrix: null columns first, then parameters of interest
    pars_of_interest = [c for c in col_full if c not in col_null]
    reordered_cols = col_null + pars_of_interest
    col_idx = [col_full.index(c) for c in reordered_cols]
    X_full = X_full[:, col_idx]
    col_full = reordered_cols

    # Indices of parameters of interest (0-based)
    index_poi = np.array(
        [i for i, c in enumerate(col_full) if c in pars_of_interest]
    )

    n_para_full = X_full.shape[1]

    # Initialize parameters
    if par_initial_full is None:
        par_initial_full = np.full(n_para_full, 0.05)
    elif len(par_initial_full) != n_para_full:
        par_initial_full = np.full(n_para_full, par_initial_full[0])

    if par_initial_null is None:
        par_initial_null = par_initial_full.copy()
    elif len(par_initial_null) != n_para_full:
        par_initial_null = np.full(n_para_full, par_initial_null[0])

    # Ensure null parameters of interest are zero
    par_initial_null[index_poi] = 0.0

    # Fit null model
    null_para = estimate_parameters_null(
        r_by_c=pic_matrix,
        design_mat=X_full,
        par_initial=par_initial_null,
        hold_zero=index_poi,
        cap_rate_vec=cap_rates,
        n_jobs=n_jobs,
    )

    # Fit full model
    full_para = estimate_parameters(
        r_by_c=pic_matrix,
        design_mat=X_full,
        par_initial=par_initial_full,
        cap_rate_vec=cap_rates,
        n_jobs=n_jobs,
    )

    # Extract convergence status and parameters
    conv_null = null_para[-1, :]
    conv_full = full_para[-1, :]
    null_coef = null_para[:-1, :]
    full_coef = full_para[:-1, :]

    # Compute p-values via gLRT
    p_values = compare_models(
        X_full=X_full,
        theta_full=full_coef,
        X_null=X_full,
        theta_null=null_coef,
        q_vec=cap_rates,
        c_by_r=pic_matrix.T,
        df_test=len(index_poi),
        n_jobs=n_jobs,
    )

    # Build feature names
    feature_names = (
        [f"f_{i + 1}" for i in range(pic_matrix.shape[0])]
    )

    p_series = pd.Series(p_values, index=feature_names, name="p_value")
    conv_df = pd.DataFrame(
        {"null": conv_null.astype(int), "full": conv_full.astype(int)},
        index=feature_names,
    )

    return PACSResult(p_values=p_series, convergence=conv_df)


def pacs_test_cumu(
    metadata: pd.DataFrame,
    formula_full: str,
    formula_null: str,
    pic_matrix: np.ndarray,
    cap_rates: np.ndarray,
    max_t: int = 2,
    n_jobs: int = 1,
    par_initial_null: np.ndarray | None = None,
    par_initial_full: np.ndarray | None = None,
) -> PACSResult:
    """Likelihood ratio test using cumulative logit model.

    Parameters
    ----------
    metadata : pd.DataFrame
        Cell metadata.
    formula_full : str
        Formula for the full model.
    formula_null : str
        Formula for the null model.
    pic_matrix : np.ndarray, shape (n_features, n_cells)
        Region-by-cell PIC count matrix (values can be 0, 1, 2, ...).
    cap_rates : np.ndarray, shape (n_cells,)
        Capturing probability for each cell.
    max_t : int
        Maximum accessibility threshold considered.
    n_jobs : int
        Number of parallel jobs.
    par_initial_null : np.ndarray, optional
        Initial null parameters.
    par_initial_full : np.ndarray, optional
        Initial full parameters.

    Returns
    -------
    PACSResult
        Test results.
    """
    # Build design matrices
    X_full_orig, col_full_orig = build_design_matrix(formula_full, metadata)
    X_null_orig, col_null_orig = build_design_matrix(formula_null, metadata)

    n_cells = X_full_orig.shape[0]

    # Stack the design matrix for cumulative logit
    # Beta part (covariates without intercept)
    X_full_beta = X_full_orig[:, 1:]  # remove intercept
    col_full_beta = col_full_orig[1:]

    # Stack beta part max_t times
    X_full_stacked = np.tile(X_full_beta, (max_t, 1))

    # Intercept part: block-diagonal identity for each threshold
    X_alpha = np.zeros((max_t * n_cells, max_t))
    for t in range(max_t):
        X_alpha[t * n_cells : (t + 1) * n_cells, t] = 1.0
    col_alpha = [f"intercept_{t + 1}" for t in range(max_t)]

    # Combined full design matrix
    X_full = np.hstack([X_alpha, X_full_stacked])
    col_full = col_alpha + col_full_beta

    # Null design matrix
    if len(col_null_orig) > 1:
        X_null_beta = X_null_orig[:, 1:]
        col_null_beta = col_null_orig[1:]
        X_null_stacked = np.tile(X_null_beta, (max_t, 1))
        X_null = np.hstack([X_alpha, X_null_stacked])
        col_null = col_alpha + col_null_beta
    else:
        X_null = X_alpha
        col_null = col_alpha

    # Reorder: null columns first, then parameters of interest
    pars_of_interest = [c for c in col_full if c not in col_null]
    reordered_cols = col_null + pars_of_interest
    col_idx = [col_full.index(c) for c in reordered_cols]
    X_full = X_full[:, col_idx]
    col_full = reordered_cols

    index_poi = np.array(
        [i for i, c in enumerate(col_full) if c in pars_of_interest]
    )

    n_para_full = X_full.shape[1]

    # Initialize parameters
    if par_initial_full is None:
        par_initial_full = np.full(n_para_full, 0.05)
    elif len(par_initial_full) != n_para_full:
        par_initial_full = np.full(n_para_full, par_initial_full[0])

    if par_initial_null is None:
        par_initial_null = par_initial_full.copy()
    elif len(par_initial_null) != n_para_full:
        par_initial_null = np.full(n_para_full, par_initial_null[0])

    par_initial_null[index_poi] = 0.0

    # Stack the PIC matrix: for each threshold t, create I(y >= t)
    pic_cumu_list = []
    for t in range(1, max_t + 1):
        pic_cumu_list.append((pic_matrix >= t).astype(np.float64))
    # Stack horizontally: each feature now has max_t * n_cells columns
    pic_stacked = np.hstack(pic_cumu_list)

    # Replicate cap_rates for stacked structure
    cap_rates_stacked = np.tile(cap_rates, max_t)

    # Fit models
    null_para = estimate_parameters_null(
        r_by_c=pic_stacked,
        design_mat=X_full,
        par_initial=par_initial_null,
        hold_zero=index_poi,
        cap_rate_vec=cap_rates_stacked,
        n_jobs=n_jobs,
    )

    full_para = estimate_parameters(
        r_by_c=pic_stacked,
        design_mat=X_full,
        par_initial=par_initial_full,
        cap_rate_vec=cap_rates_stacked,
        n_jobs=n_jobs,
    )

    # Extract convergence and coefficients
    conv_null = null_para[-1, :]
    conv_full = full_para[-1, :]
    null_coef = null_para[:-1, :]
    full_coef = full_para[:-1, :]

    # P-values
    p_values = compare_models(
        X_full=X_full,
        theta_full=full_coef,
        X_null=X_full,
        theta_null=null_coef,
        q_vec=cap_rates_stacked,
        c_by_r=pic_stacked.T,
        df_test=len(index_poi),
        n_jobs=n_jobs,
    )

    feature_names = [f"f_{i + 1}" for i in range(pic_stacked.shape[0])]

    p_series = pd.Series(p_values, index=feature_names, name="p_value")
    conv_df = pd.DataFrame(
        {"null": conv_null.astype(int), "full": conv_full.astype(int)},
        index=feature_names,
    )

    return PACSResult(p_values=p_series, convergence=conv_df)


def pacs_test(
    pic_matrix: np.ndarray | sp.spmatrix,
    metadata: pd.DataFrame,
    formula_full: str,
    formula_null: str,
    cap_rates: np.ndarray,
    *,
    t_proportion_cutoff: float = 0.25,
    n_peaks_per_round: int | None = None,
    n_jobs: int = 1,
    verbose: bool = True,
    par_initial_null: np.ndarray | None = None,
    par_initial_full: np.ndarray | None = None,
) -> PACSResult:
    """PACS test with automatic model selection.

    Automatically selects between binary logit and cumulative logit models
    for each feature based on the proportion of counts >= 2.

    Parameters
    ----------
    pic_matrix : np.ndarray or scipy.sparse matrix, shape (n_features, n_cells)
        Region-by-cell PIC count matrix.
    metadata : pd.DataFrame
        Cell metadata with columns referenced in the formulas.
    formula_full : str
        R-style formula for the full model, e.g., ``'~ cell_type + batch'``.
    formula_null : str
        R-style formula for the null model, e.g., ``'~ batch'``.
    cap_rates : np.ndarray, shape (n_cells,)
        Capturing probability for each cell.
    t_proportion_cutoff : float
        Threshold for selecting cumulative vs. logit model. Features with
        proportion of counts >= 2 above this use cumulative logit.
    n_peaks_per_round : int, optional
        Number of peaks per batch. Auto-determined if None.
    n_jobs : int
        Number of parallel jobs.
    verbose : bool
        Whether to print progress.
    par_initial_null : np.ndarray, optional
        Initial null parameters.
    par_initial_full : np.ndarray, optional
        Initial full parameters.

    Returns
    -------
    PACSResult
        Test results with p-values and convergence status for all features.
    """
    is_sparse = sp.issparse(pic_matrix)

    n_features, n_cells = pic_matrix.shape

    # Validate inputs
    if metadata.shape[0] != n_cells:
        raise ValueError(
            f"Number of cells in metadata ({metadata.shape[0]}) does not "
            f"match data matrix ({n_cells})"
        )
    if len(cap_rates) != n_cells:
        raise ValueError(
            f"Length of cap_rates ({len(cap_rates)}) does not match "
            f"data matrix ({n_cells})"
        )

    # Generate feature names if needed
    feature_names = np.array([f"f_{i + 1}" for i in range(n_features)])

    # Determine batch size
    if n_peaks_per_round is None:
        n_peaks_per_round = min(2**29 // n_cells, n_features)

    # Compute proportion of counts >= 2 per feature
    if is_sparse:
        pic_sp = sp.csr_matrix(pic_matrix)

        # Binarized matrix (any nonzero -> 1)
        pic_bin = pic_sp.copy()
        pic_bin.data = np.ones_like(pic_bin.data)

        # Counts >= 2 matrix
        pic_2 = pic_sp.copy()
        pic_2.data[pic_2.data == 1] = 0
        pic_2.eliminate_zeros()
        pic_2.data = np.ones_like(pic_2.data)

        rs = np.asarray(pic_bin.sum(axis=1)).ravel()
        rs2 = np.asarray(pic_2.sum(axis=1)).ravel()
    else:
        pic_bin = (pic_matrix > 0).astype(np.float64)
        pic_2 = (pic_matrix >= 2).astype(np.float64)

        rs = pic_bin.sum(axis=1)
        rs2 = pic_2.sum(axis=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        p_2 = np.where(rs > 0, rs2 / rs, 0.0)

    # Split features
    cumu_mask = p_2 >= t_proportion_cutoff
    logit_mask = ~cumu_mask
    n_cumu = int(np.sum(cumu_mask))
    n_logit = int(np.sum(logit_mask))

    if verbose:
        print(f"{n_cumu} peaks use cumulative logit model")
        print(f"{n_logit} peaks use binary logit model")

    # Initialize result arrays
    all_p_values = np.full(n_features, np.nan)
    all_conv_null = np.full(n_features, np.nan)
    all_conv_full = np.full(n_features, np.nan)

    # --- Cumulative logit part ---
    if n_cumu > 0:
        cumu_indices = np.where(cumu_mask)[0]
        n_iters = int(np.ceil(n_cumu / n_peaks_per_round))

        for jj in range(n_iters):
            start = jj * n_peaks_per_round
            end = min(n_cumu, (jj + 1) * n_peaks_per_round)
            batch_idx = cumu_indices[start:end]

            if is_sparse:
                pic_dense = np.asarray(pic_matrix[batch_idx, :].todense())
            else:
                pic_dense = pic_matrix[batch_idx, :]

            result = pacs_test_cumu(
                metadata=metadata,
                formula_full=formula_full,
                formula_null=formula_null,
                pic_matrix=pic_dense,
                cap_rates=cap_rates,
                n_jobs=n_jobs,
                par_initial_null=par_initial_null,
                par_initial_full=par_initial_full,
            )

            all_p_values[batch_idx] = result.p_values.values
            all_conv_null[batch_idx] = result.convergence["null"].values
            all_conv_full[batch_idx] = result.convergence["full"].values

    # --- Binary logit part ---
    if n_logit > 0:
        logit_indices = np.where(logit_mask)[0]
        n_iters_b = int(np.ceil(n_logit / n_peaks_per_round))

        for jj in range(n_iters_b):
            start = jj * n_peaks_per_round
            end = min(n_logit, (jj + 1) * n_peaks_per_round)
            batch_idx = logit_indices[start:end]

            if is_sparse:
                pic_dense = np.asarray(pic_matrix[batch_idx, :].todense())
            else:
                pic_dense = pic_matrix[batch_idx, :]

            # Binarize for logit model
            pic_dense = (pic_dense > 0).astype(np.float64)

            result = pacs_test_logit(
                metadata=metadata,
                formula_full=formula_full,
                formula_null=formula_null,
                pic_matrix=pic_dense,
                cap_rates=cap_rates,
                n_jobs=n_jobs,
                par_initial_null=par_initial_null,
                par_initial_full=par_initial_full,
            )

            all_p_values[batch_idx] = result.p_values.values
            all_conv_null[batch_idx] = result.convergence["null"].values
            all_conv_full[batch_idx] = result.convergence["full"].values

    p_series = pd.Series(all_p_values, index=feature_names, name="p_value")
    conv_df = pd.DataFrame(
        {"null": all_conv_null.astype(int), "full": all_conv_full.astype(int)},
        index=feature_names,
    )

    return PACSResult(p_values=p_series, convergence=conv_df)
