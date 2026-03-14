"""PACS: Probabilistic model for Accessible Chromatin data in Single Cell.

A Python toolkit for snATAC-seq data analysis, providing differential
accessible region identification and cell type annotation.

Reference
---------
Miao et al. "Model-based compound hypothesis testing for snATAC-seq data
with PACS" bioRxiv (2023).
"""

__version__ = "0.1.0"

from ._types import AnnotationResult, ConvergenceStatus, PACSResult
from .annotation import estimate_label, estimate_label_selected_peaks
from .cct import cauchy_combination_test
from .testing import pacs_test, pacs_test_cumu, pacs_test_logit

__all__ = [
    "pacs_test",
    "pacs_test_logit",
    "pacs_test_cumu",
    "estimate_label",
    "estimate_label_selected_peaks",
    "cauchy_combination_test",
    "PACSResult",
    "AnnotationResult",
    "ConvergenceStatus",
]
