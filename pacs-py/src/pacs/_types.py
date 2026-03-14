"""Result types and enums for PACS."""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum

import numpy as np
import pandas as pd


class ConvergenceStatus(IntEnum):
    """Convergence status codes for parameter estimation."""

    CONVERGED = 1
    SINGULAR_MATRIX = 2
    MAX_ITER_REACHED = 3
    LOSS_NOT_INCREASING = 4


@dataclass
class PACSResult:
    """Result of a PACS differential accessibility test.

    Attributes
    ----------
    p_values : pd.Series
        P-values for each feature (region/peak), indexed by feature name.
    convergence : pd.DataFrame
        Convergence status with columns 'null' and 'full', indexed by feature.
        Values correspond to ConvergenceStatus codes.
    """

    p_values: pd.Series
    convergence: pd.DataFrame


@dataclass
class AnnotationResult:
    """Result of cell type annotation.

    Attributes
    ----------
    log_likelihoods : pd.DataFrame
        Log-likelihood matrix (cells x cell_types).
    predicted_labels : pd.Series
        Most likely cell type for each cell.
    """

    log_likelihoods: pd.DataFrame
    predicted_labels: pd.Series
