"""Utility functions for PACS."""

from __future__ import annotations

import numpy as np
import pandas as pd
from formulaic import model_matrix


def build_design_matrix(
    formula: str, data: pd.DataFrame
) -> tuple[np.ndarray, list[str]]:
    """Convert an R-style formula string to a design matrix.

    Parameters
    ----------
    formula : str
        R-style formula, e.g., ``'~ cell_type + batch'``.
    data : pd.DataFrame
        Metadata with columns referenced in the formula.

    Returns
    -------
    matrix : np.ndarray
        Design matrix with intercept and dummy variables.
    col_names : list[str]
        Column names of the design matrix.
    """
    dm = model_matrix(formula, data)
    return np.asarray(dm, dtype=np.float64), list(dm.columns)
