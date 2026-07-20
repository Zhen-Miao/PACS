# Cauchy Combination Test accelerated

This is adopted from the Liu and Xie JASA 2020 publication and their
GitHub page. Please also cite their initial publication if you used this
part of codes. Thaks!

## Usage

``` r
CCT_internal_horizontal(pval_mat)
```

## Arguments

- pval_mat:

  A matrix of p values where each row is a feature and each column is a
  p value from one individual test, so that we want to compute
  aggregated p values for each row

## Value

A vector of p value for each feature after Cauchy combination
