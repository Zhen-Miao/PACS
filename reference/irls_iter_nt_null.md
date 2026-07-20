# Iteratively reweighted least squares function for null – newton method for the null model

Iteratively reweighted least squares function for null – newton method
for the null model

## Usage

``` r
irls_iter_nt_null(
  y_vec,
  xdumm,
  theta_estimated,
  hold_zero,
  tolerance = .Machine$double.eps^0.5,
  stop_criteria = 1e-05,
  q_vec
)
```

## Arguments

- y_vec:

  Observed counts

- xdumm:

  Design matrix

- theta_estimated:

  Estimated coefficient values

- hold_zero:

  Which element to hold zero in the null model

- tolerance:

  Numeric tolerance for matrix to be considered not revertable, default
  = .Machine\$double.eps^0.5

- stop_criteria:

  Resolution to stop, default = 1e-6

- q_vec:

  Capturing probability vector

## Value

Converged coefficient estimate for the null model
