# Use gLRT to compare two models and obtain p values

Use gLRT to compare two models and obtain p values

## Usage

``` r
compare_models(
  x_full,
  theta_estimated_full,
  x_null,
  theta_estimated_null,
  q_vec,
  c_by_r,
  df_test,
  mc.cores = 1
)
```

## Arguments

- x_full:

  Design matrix for the full model

- theta_estimated_full:

  Estimated coefficients for the full model

- x_null:

  Design matrix for the null model

- theta_estimated_null:

  Estimated coefficients for the null model

- q_vec:

  A vector of capturing probability for each cell

- c_by_r:

  Cell by region matrix

- df_test:

  Degree of freedom of the test

- mc.cores:

  Number of cores in multi-core computing

## Value

A vector of p values based on gLRT
