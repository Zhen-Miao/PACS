## T=2 parameter recovery: simulate from known (alpha, beta) under Option B,
## fit, and assert MLE within tolerance of truth across replicates.
##
## Soft assertions (expect_lt) are used so flakiness is bounded; thresholds
## are sized for the n=400, 30-replicate setup below.

test_that("MLE recovers alpha and beta at T=2", {
  skip_on_cran()
  n_rep <- 30L
  n <- 400L
  alpha <- c(0.8, -0.4)
  beta <- c(0.6, -0.3)

  alpha_hat <- matrix(NA_real_, nrow = n_rep, ncol = 2L)
  beta_hat <- matrix(NA_real_, nrow = n_rep, ncol = 2L)
  conv <- integer(n_rep)

  for (r in seq_len(n_rep)) {
    fx <- make_cumu_fixture(n = n, p = 2L, T = 2L,
                            alpha = alpha, beta = beta, seed = 100L + r)
    M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                            q = fx$q, capture = "B", seed = 200L + r)
    theta_init <- PACS:::warm_start_theta(M = M, q = fx$q, T = 2L, p_beta = 2L)
    res <- PACS:::irls_iter_cumu(M_vec = M, X = fx$X,
                                 theta_estimated = theta_init,
                                 q_vec = fx$q, T = 2L)
    conv[r] <- res[length(res)]
    th <- res[seq_len(length(res) - 1L)]
    if (conv[r] == 1L) {
      ## Map atilde back to alpha for comparison.
      alpha_hat[r, ] <- PACS:::alpha_from_atilde(th[1:2])
      beta_hat[r, ] <- th[3:4]
    }
  }

  ok <- conv == 1L
  expect_gt(mean(ok), 0.85)  ## convergence rate

  bias_alpha <- colMeans(alpha_hat[ok, , drop = FALSE]) - alpha
  bias_beta <- colMeans(beta_hat[ok, , drop = FALSE]) - beta

  ## Mean absolute bias should be small with n=400 and 30 replicates.
  expect_lt(max(abs(bias_alpha)), 0.15)
  expect_lt(max(abs(bias_beta)), 0.15)
})
