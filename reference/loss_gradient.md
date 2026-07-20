# Gradient of loss function without Firth prior

Gradient of loss function without Firth prior

## Usage

``` r
loss_gradient(xdumm, p_bg, q_vec, y_vec)
```

## Arguments

- xdumm:

  A design matrix, X with dummy variables

- p_bg:

  A value of open probability

- q_vec:

  A vector of capturing rates

- y_vec:

  A vector of observed open (1) or close (0) state for each cell, must
  match the length of q

## Value

The first order derivative of the loss function without Firth prior, a
vector
