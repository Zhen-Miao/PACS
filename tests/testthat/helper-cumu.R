## Shared fixtures for cumulative-logit tests.

## Default seed-controlled fixture: n cells, p covariates, T thresholds,
## known (alpha, beta), capture rates uniformly in [q_lo, q_hi].
make_cumu_fixture <- function(n = 400L, p = 2L, T = 2L,
                              alpha = c(0.8, -0.4),
                              beta = c(0.6, -0.3),
                              q_lo = 0.4, q_hi = 0.9, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  q <- runif(n, q_lo, q_hi)
  list(X = X, q = q, alpha = alpha[seq_len(T)], beta = beta[seq_len(p)],
       T = T, n = n, p = p)
}

## Inverse of PACS:::alpha_from_atilde, written independently of it so that
## tests parameterise by the ordered thresholds rather than by the internal
## unconstrained coordinates.
atilde_from_alpha_test <- function(alpha) {
  if (length(alpha) == 1L) return(alpha)
  c(alpha[1L], log(alpha[-length(alpha)] - alpha[-1L]))
}

## Build (alpha_1, beta) -> theta-binary equivalent vector for parity tests.
binary_theta_from_cumu <- function(theta_cumu, T) {
  ## In the cumu code at T=1, theta = (atilde_1 = alpha_1, beta...).
  ## In binary code, theta = (intercept, beta...) with X having an intercept
  ## column. So they map directly: binary_theta = theta_cumu.
  stopifnot(T == 1L)
  theta_cumu
}
