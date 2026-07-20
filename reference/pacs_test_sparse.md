# Profile likelihood ratio test with PACS

Profile likelihood ratio test with PACS

## Usage

``` r
pacs_test_sparse(
  covariate_meta.data,
  formula_full,
  formula_null,
  pic_matrix,
  n_peaks_per_round = NULL,
  T_proportion_cutoff = 0.2,
  cap_rates,
  par_initial_null = NULL,
  par_initial_full = NULL,
  n_cores = 1,
  verbose = TRUE
)
```

## Arguments

- covariate_meta.data:

  A data.frame with columns representing the covariates and rows
  representing cells

- formula_full:

  A formula object representing the full model. For example, ~
  cell_type + batch

- formula_null:

  A formula object representing the null model. For example, ~ batch

- pic_matrix:

  The input region-by-cell PIC matrix

- n_peaks_per_round:

  The number of peaks tested in each round. Default = NULL, which will
  be determined by the maximum allowable size of matrix in R.

- T_proportion_cutoff:

  To determine the maximum value of T, we set a criteria based on the
  proportion of reads that are \>= T. This parameter set this
  proportion. Default = 0.2

- cap_rates:

  A vector of capturing probability for each cell

- par_initial_null:

  Initialized values of estimated parameters for the null model, we do
  not need to specify unless there are reasons to do so. Default = NULL

- par_initial_full:

  Initialized values of estimated parameters for the null model, we do
  not need to specify unless there are reasons to do so. Default = NULL

- n_cores:

  number of cores for multi-core computation

- verbose:

  Whether to print progress, default = TRUE

## Value

A list of two elements, pacs_converged is a vector of length 2\*n_peaks
representing the convergence status of the peak in the null and full
model and pacs_p_val is a vector of length n_peaks representing the p
values for each peak.
