"""Cell type annotation via likelihood-based classification.

Ported from ``cell_type_annotation.R``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ._types import AnnotationResult


def estimate_label(
    r_by_t: np.ndarray,
    in_r_by_c: np.ndarray | sp.spmatrix,
    *,
    alpha: float = 1.0,
) -> AnnotationResult:
    """Estimate cell type labels for unannotated cells.

    Uses a reference region-by-cell-type probability matrix to classify
    cells based on their accessibility profile.

    Parameters
    ----------
    r_by_t : np.ndarray, shape (n_regions, n_cell_types)
        Region-by-cell-type probability matrix from an annotated dataset.
        Column names should be cell type labels.
    in_r_by_c : np.ndarray or scipy.sparse matrix, shape (n_regions, n_cells)
        Binary region-by-cell matrix (unannotated). Values must be 0 or 1.
    alpha : float
        Weight for negative peaks (regions where cell is not accessible).

    Returns
    -------
    AnnotationResult
        Log-likelihoods and predicted cell type labels.
    """
    # Validate binary matrix
    if sp.issparse(in_r_by_c):
        data = in_r_by_c.data
        if np.any((data != 0) & (data != 1)):
            raise ValueError("Input matrix must be binary (0/1).")
        in_r_by_c_dense = np.asarray(in_r_by_c.todense())
    else:
        if np.any((in_r_by_c != 0) & (in_r_by_c != 1)):
            raise ValueError("Input matrix must be binary (0/1).")
        in_r_by_c_dense = np.asarray(in_r_by_c, dtype=np.float64)

    r_by_t = np.asarray(r_by_t, dtype=np.float64)

    n_regions, n_cells = in_r_by_c_dense.shape
    n_cell_types = r_by_t.shape[1]

    # Estimate capturing rates per cell per cell type
    n_reads = in_r_by_c_dense.sum(axis=0)  # (n_cells,)
    sum_prob = r_by_t.sum(axis=0)  # (n_cell_types,)

    # cap_rate[ct, cell] = n_reads[cell] / sum_prob[ct]
    cap_rate_mat = n_reads[np.newaxis, :] / sum_prob[:, np.newaxis]
    cap_rate_mat = np.clip(cap_rate_mat, 0.00005, 0.9995)

    # Process in batches if matrix is too large
    max_elements = 2**31 - 1
    total_elements = n_regions * n_cells
    if total_elements > max_elements:
        split_by = int(np.ceil(total_elements / max_elements))
        batch_size = int(np.ceil(n_cells / split_by))
    else:
        split_by = 1
        batch_size = n_cells

    all_results = []

    for batch_start in range(0, n_cells, batch_size):
        batch_end = min(batch_start + batch_size, n_cells)
        x_batch = in_r_by_c_dense[:, batch_start:batch_end]
        n_batch = batch_end - batch_start

        # Compute log-likelihoods for each cell in batch
        batch_ll = np.zeros((n_batch, n_cell_types))

        for ct in range(n_cell_types):
            q_ct = cap_rate_mat[ct, batch_start:batch_end]  # (n_batch,)
            # pq[r, c] = r_by_t[r, ct] * q_ct[c]
            pq = r_by_t[:, ct : ct + 1] * q_ct[np.newaxis, :]
            log_pq = np.log(pq)
            log_1_pq = np.log1p(-pq)

            # ll[c] = sum_r(x[r,c] * log(pq[r,c]) + (1-x[r,c]) * log(1-pq[r,c]) * alpha)
            batch_ll[:, ct] = np.sum(
                x_batch * log_pq + (1.0 - x_batch) * log_1_pq * alpha,
                axis=0,
            )

        all_results.append(batch_ll)

    esti_mat = np.vstack(all_results)

    # Create result with proper labels
    cell_names = [f"cell_{i + 1}" for i in range(n_cells)]
    ct_names = [f"ct_{i + 1}" for i in range(n_cell_types)]

    ll_df = pd.DataFrame(esti_mat, index=cell_names, columns=ct_names)
    predicted = pd.Series(
        ll_df.columns[np.argmax(esti_mat, axis=1)],
        index=cell_names,
        name="predicted_label",
    )

    return AnnotationResult(log_likelihoods=ll_df, predicted_labels=predicted)


def estimate_label_selected_peaks(
    r_by_t: np.ndarray,
    in_r_by_c: np.ndarray | sp.spmatrix,
    peaks_sel: np.ndarray | list,
    *,
    alpha: float = 1.0,
) -> AnnotationResult:
    """Estimate cell type labels using selected informative peaks.

    Parameters
    ----------
    r_by_t : np.ndarray, shape (n_regions, n_cell_types)
        Region-by-cell-type probability matrix.
    in_r_by_c : np.ndarray or scipy.sparse matrix, shape (n_regions, n_cells)
        Binary region-by-cell matrix.
    peaks_sel : array-like
        Indices of selected informative peaks.
    alpha : float
        Weight for negative peaks.

    Returns
    -------
    AnnotationResult
        Log-likelihoods and predicted labels.
    """
    peaks_sel = np.asarray(peaks_sel)
    return estimate_label(
        r_by_t=r_by_t[peaks_sel, :],
        in_r_by_c=in_r_by_c[peaks_sel, :],
        alpha=alpha,
    )
