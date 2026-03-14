# PACS-py: Python Implementation of PACS

**PACS** (Probabilistic model for Accessible Chromatin data in Single Cell) is a toolkit for snATAC-seq data analysis. This is the Python implementation, designed to integrate with the Python single-cell ecosystem (AnnData, scanpy).

## Reference

Miao et al. "Model-based compound hypothesis testing for snATAC-seq data with PACS" *bioRxiv* (2023).

## Installation

```bash
pip install -e .
```

With AnnData integration:

```bash
pip install -e ".[anndata]"
```

For development:

```bash
pip install -e ".[dev]"
```

## Quick Start

### Differential Accessible Region (DAR) Identification

```python
import numpy as np
import pandas as pd
from pacs import pacs_test

# pic_matrix: regions x cells PIC count matrix
# metadata: DataFrame with cell covariates
# cap_rates: capturing probability per cell

result = pacs_test(
    pic_matrix=pic_matrix,
    metadata=metadata,
    formula_full="~ cell_type + batch",
    formula_null="~ batch",
    cap_rates=cap_rates,
    n_jobs=4,
)

# Access results
significant = result.p_values[result.p_values < 0.05]
print(f"Found {len(significant)} DARs")
```

### Cell Type Annotation

```python
from pacs import estimate_label

# r_by_t: regions x cell_types probability matrix (from reference)
# pic_binary: binary regions x cells matrix (query data)

result = estimate_label(r_by_t, pic_binary)
print(result.predicted_labels)
```

### AnnData Integration

```python
from pacs.io import pacs_test_anndata

# adata.obs must contain covariates and a 'cap_rates' column
result = pacs_test_anndata(
    adata,
    formula_full="~ cell_type",
    formula_null="~ 1",
    n_jobs=4,
)
# Results stored in adata.var['pacs_pvalue']
```

### Cauchy Combination Test

```python
from pacs import cauchy_combination_test

# Combine p-values from multiple tests (features x tests)
combined_pvals = cauchy_combination_test(p_value_matrix)
```

## API Reference

### Main Functions

| Function | Description |
|----------|-------------|
| `pacs_test()` | Auto-selects logit or cumulative logit model per feature |
| `pacs_test_logit()` | Binary logit model for DAR testing |
| `pacs_test_cumu()` | Cumulative logit model for count data |
| `estimate_label()` | Cell type annotation from reference |
| `estimate_label_selected_peaks()` | Annotation using selected peaks |
| `cauchy_combination_test()` | P-value aggregation across tests |

### Result Types

- **`PACSResult`**: Contains `p_values` (pd.Series) and `convergence` (pd.DataFrame)
- **`AnnotationResult`**: Contains `log_likelihoods` (pd.DataFrame) and `predicted_labels` (pd.Series)
- **`ConvergenceStatus`**: Enum with `CONVERGED`, `SINGULAR_MATRIX`, `MAX_ITER_REACHED`, `LOSS_NOT_INCREASING`

## Dependencies

- numpy >= 1.21
- scipy >= 1.7
- pandas >= 1.3
- formulaic >= 0.6
- joblib >= 1.1
- anndata >= 0.8 (optional, for AnnData integration)
