"""AnnData integration for PACS.

Provides convenience wrappers that accept AnnData objects as input,
extracting the PIC matrix, metadata, and capturing rates automatically.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ._types import AnnotationResult, PACSResult
from .annotation import estimate_label
from .testing import pacs_test

if TYPE_CHECKING:
    import anndata


def pacs_test_anndata(
    adata: anndata.AnnData,
    formula_full: str,
    formula_null: str,
    *,
    cap_rates_key: str = "cap_rates",
    layer: str | None = None,
    n_jobs: int = 1,
    t_proportion_cutoff: float = 0.25,
    verbose: bool = True,
    store_result: bool = True,
) -> PACSResult:
    """Run PACS differential test on an AnnData object.

    AnnData stores data as cells-by-features (obs x var). This function
    transposes to features-by-cells as required by PACS.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix. ``adata.obs`` must contain the covariates
        referenced in the formulas and a column with capturing rates.
    formula_full : str
        R-style formula for the full model.
    formula_null : str
        R-style formula for the null model.
    cap_rates_key : str
        Column name in ``adata.obs`` for capturing rates.
    layer : str, optional
        Layer to use. If None, uses ``adata.X``.
    n_jobs : int
        Number of parallel jobs.
    t_proportion_cutoff : float
        Threshold for model selection.
    verbose : bool
        Whether to print progress.
    store_result : bool
        If True, store p-values in ``adata.var['pacs_pvalue']``.

    Returns
    -------
    PACSResult
        Test results.
    """
    # Extract matrix (cells x features -> transpose to features x cells)
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    # Transpose: AnnData is cells x features, PACS needs features x cells
    pic_matrix = X.T

    # Metadata
    metadata = adata.obs.copy()

    # Capturing rates
    if cap_rates_key not in metadata.columns:
        raise ValueError(
            f"Capturing rates column '{cap_rates_key}' not found in adata.obs. "
            f"Available columns: {list(metadata.columns)}"
        )
    cap_rates = np.asarray(metadata[cap_rates_key], dtype=np.float64)

    result = pacs_test(
        pic_matrix=pic_matrix,
        metadata=metadata,
        formula_full=formula_full,
        formula_null=formula_null,
        cap_rates=cap_rates,
        t_proportion_cutoff=t_proportion_cutoff,
        n_jobs=n_jobs,
        verbose=verbose,
    )

    # Update feature names to match adata.var_names
    if adata.var_names is not None and len(adata.var_names) > 0:
        result.p_values.index = adata.var_names
        result.convergence.index = adata.var_names

    # Store results back
    if store_result:
        adata.var["pacs_pvalue"] = result.p_values.values
        adata.var["pacs_conv_null"] = result.convergence["null"].values
        adata.var["pacs_conv_full"] = result.convergence["full"].values

    return result


def annotate_anndata(
    adata: anndata.AnnData,
    r_by_t: np.ndarray,
    *,
    layer: str | None = None,
    alpha: float = 1.0,
    store_result: bool = True,
) -> AnnotationResult:
    """Run cell type annotation on an AnnData object.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix. Must contain a binary accessibility matrix.
    r_by_t : np.ndarray, shape (n_regions, n_cell_types)
        Region-by-cell-type probability matrix from a reference dataset.
    layer : str, optional
        Layer to use. If None, uses ``adata.X``.
    alpha : float
        Weight for negative peaks.
    store_result : bool
        If True, store predicted labels in ``adata.obs['pacs_label']``.

    Returns
    -------
    AnnotationResult
        Annotation results.
    """
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    # Transpose: cells x features -> features x cells
    pic_matrix = X.T

    result = estimate_label(r_by_t=r_by_t, in_r_by_c=pic_matrix, alpha=alpha)

    # Update cell names
    if adata.obs_names is not None and len(adata.obs_names) > 0:
        result.log_likelihoods.index = adata.obs_names
        result.predicted_labels.index = adata.obs_names

    if store_result:
        adata.obs["pacs_label"] = result.predicted_labels.values

    return result
