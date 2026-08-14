# likelihood ratio test with PACS – cumulative logit

likelihood ratio test with PACS – cumulative logit

## Usage

``` r
pacs_test_cumu(
  covariate_meta.data,
  formula_full,
  formula_null,
  pic_matrix,
  max_T = 2,
  cap_rates,
  par_initial_null = NULL,
  par_initial_full = NULL,
  n_cores = 1,
  method = c("stacked", "exact")
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

- max_T:

  The maximum accessibility category considered, default = 2. For
  `method = "exact"`, observed counts greater than `max_T` are
  top-coded, so the final category means `max_T` or more.

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

- method:

  One of `"stacked"` (default, back-compat) or `"exact"`. `"stacked"`
  uses the original stack-and-treat-as-binary approximation. `"exact"`
  uses the proper cumulative-logit likelihood with capture-rate
  correction (Option B in `notes/cumulative_logit_math.md`) and an
  unpenalized maximum-likelihood fit. Exact-path p-values are ordinary
  likelihood-ratio tests and are `NA` for non-converged or boundary
  fits.

## Value

A list of two elements. `pacs_converged` has length `2 * n_peaks`, with
null-fit statuses followed by full-fit statuses. For the exact path,
status 1 means converged, 2 means a singular or non-finite scoring
system, 3 means the iteration limit was reached, 4 means step-halving
found no acceptable update, and 5 means the supplied starting value had
a non-finite log-likelihood. `pacs_p_val` contains one p-value per peak;
exact-path inference is `NA` unless both statuses are 1 and both fits
are interior.

## Details

The unpenalized exact MLE can fail to exist or have singular information
for very sparse peaks. In review simulations, the exact-path `NA` rate
rose from 0% for dense peaks to 3% for moderately sparse peaks
(`n = 300`, `alpha = c(-2.5, -4)`), 38% for still sparser peaks
(`n = 300`, `alpha = c(-3.5, -5)`), and 46% with fewer cells (`n = 100`,
`alpha = c(-2.5, -4)`). These rates are scenario-specific, not general
guarantees. Withholding a p-value is a conservative failure policy for
the affected peak—it avoids turning a failed fit into a false
positive—but the resulting loss of analyzable peaks reduces power.
