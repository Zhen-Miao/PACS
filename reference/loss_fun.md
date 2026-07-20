# Calculate Loss Function without Firth prior

Calculate Loss Function without Firth prior

## Usage

``` r
loss_fun(p_bg, q_vec, y_vec)
```

## Arguments

- p_bg:

  A vector of open probability for cell type group g

- q_vec:

  A vector of capturing rates

- y_vec:

  A vector of observed open (1) or close (0) state for each cell, must
  match the length of q

## Value

Log loss value without Firth prior
