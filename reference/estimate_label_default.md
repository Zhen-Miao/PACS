# Estimate cell type label for new data

Estimate cell type label for new data

## Usage

``` r
estimate_label_default(r_by_t, in_r_by_c, alpha = 1)
```

## Arguments

- r_by_t:

  Region-by cell type matrix generated from a annotated dataset

- in_r_by_c:

  Input (unannotated) region by cell matrix

- alpha:

  Weight for negative peaks, default = 1

## Value

A matrix of cell types by cell matrix, with elements representing
probability of being in that cell type
