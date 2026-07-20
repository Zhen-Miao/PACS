# Estimate parameters for the null model

Estimate parameters for the null model

## Usage

``` r
estimate_parameters_null(
  r_by_c,
  design_mat,
  par_initial,
  hold_zero,
  cap_rate_vec,
  mc.cores = 1
)
```

## Arguments

- r_by_c:

  Region by cell matrix

- design_mat:

  Desing matrix

- par_initial:

  Initialized parameters

- hold_zero:

  Which parameter to hold zero during estimation

- cap_rate_vec:

  A vector of capturing probability

- mc.cores:

  Number of multi-cores, default = 1

## Value

Estimated coefficient for the null model
