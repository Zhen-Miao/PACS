## Simulator for cumulative-logit PACS data with capture rate.
##
## Internal helper used by tests/testthat/ and (optionally) vignettes.
## Not exported.


## Sample M_i from the Option B observed-count distribution.
##
## With probability q_i, draw Y_i from the proportional-odds multinomial
## defined by (alpha, beta), and set M_i = Y_i.
## With probability 1 - q_i, set M_i = 0.
##
## Args:
##   X      : n × p covariate matrix (no intercept; alpha plays that role).
##   alpha  : length-T threshold vector with alpha_1 >= ... >= alpha_T.
##   beta   : length-p slope vector.
##   q      : length-n capture rates.
##   capture: "B" only (Option A is deferred per the math note's plan).
##
## Returns: integer vector of length n.
#' @importFrom stats runif
#' @noRd
simulate_cumu_pacs <- function(X, alpha, beta, q, capture = "B", seed = NULL) {
  if (capture != "B") {
    stop("Only capture = 'B' is implemented; Option A (thinning) is deferred.")
  }
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(X)
  p <- ncol(X)
  T <- length(alpha)
  if (length(beta) != p) {
    stop("length(beta) must equal ncol(X).")
  }
  if (any(diff(alpha) > 0)) {
    stop("alpha must satisfy alpha_1 >= alpha_2 >= ... >= alpha_T.")
  }
  if (length(q) != n) {
    stop("length(q) must equal nrow(X).")
  }

  xb <- as.numeric(X %*% beta)
  ## Cumulative probabilities Pr(Y >= t | x) for t = 1..T.
  P_ge <- matrix(0, nrow = n, ncol = T)
  for (t in 1:T) {
    P_ge[, t] <- 1 / (1 + exp(-(alpha[t] + xb)))
  }
  ## Category probabilities pi_{ik} for k = 0..T.
  pi_mat <- matrix(0, nrow = n, ncol = T + 1L)
  pi_mat[, 1] <- 1 - P_ge[, 1]                       ## k = 0
  if (T >= 2L) {
    for (k in 1:(T - 1L)) {
      pi_mat[, k + 1L] <- P_ge[, k] - P_ge[, k + 1L] ## k = 1..T-1
    }
  }
  pi_mat[, T + 1L] <- P_ge[, T]                      ## k = T

  ## Sample Y_i from the per-row categorical, then drop to 0 with prob 1-q_i.
  Y <- integer(n)
  for (i in seq_len(n)) {
    Y[i] <- sample.int(T + 1L, size = 1L, prob = pi_mat[i, ]) - 1L
  }
  drop <- runif(n) > q
  M <- Y
  M[drop] <- 0L
  M
}
