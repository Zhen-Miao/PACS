# Estimate parameters for the full model

Estimate parameters for the full model

## Usage

``` r
estimate_parameters(
  r_by_c,
  design_mat,
  par_initial,
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

- cap_rate_vec:

  A vector of capturing probability

- mc.cores:

  Number of multi-cores, default = 1

## Value

Estimated coefficient for the full model
