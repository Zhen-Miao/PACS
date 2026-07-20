# Estimate cell type label for new data using all peaks as input

We will identify relevant peaks

## Usage

``` r
estimate_label_no_cap_rate_all_pks(r_by_t, in_r_by_c, pks_sel, alpha = 1)
```

## Arguments

- r_by_t:

  Region-by cell type matrix generated from a annotated dataset

- in_r_by_c:

  Input (unannotated) region by cell matrix

- pks_sel:

  selected peaks as informative for cell type label prediction

- alpha:

  Weight for negative peaks, default = 1

## Value

A matrix of cell types by cell matrix, with elements representing
probability of being in that cell type
