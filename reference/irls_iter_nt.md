# Iteratively reweighted least squares function – newton methods

Iteratively reweighted least squares function – newton methods

## Usage

``` r
irls_iter_nt(
  y_vec,
  xdumm,
  theta_estimated,
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

- tolerance:

  Numeric tolerance for matrix to be considered not revertable, default
  = .Machine\$double.eps^0.5

- stop_criteria:

  Resolution to stop, default = 1e-6

- q_vec:

  Capturing probability vector

## Value

Converged coefficient estimate
