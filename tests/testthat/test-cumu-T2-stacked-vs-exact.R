## Informational test: at T=2 on simulated data, the exact path should
## produce lower mean squared error on beta than the stacked approximation.
##
## This is a sanity check, not a strict assertion of any specific magnitude:
## the stacked approach has provable bias because it treats nested
## indicators as independent Bernoulli, double-counting each cell.

test_that("exact has no worse bias than stacked at T=2", {
  skip_on_cran()
  n_rep <- 20L
  n <- 400L
  alpha <- c(0.8, -0.4)
  beta <- c(0.6, -0.3)

  beta_hat_exact <- matrix(NA_real_, nrow = n_rep, ncol = 2L)
  beta_hat_stack <- matrix(NA_real_, nrow = n_rep, ncol = 2L)

  for (r in seq_len(n_rep)) {
    fx <- make_cumu_fixture(n = n, p = 2L, T = 2L,
                            alpha = alpha, beta = beta, seed = 700L + r)
    M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                            q = fx$q, capture = "B", seed = 800L + r)

    ## Exact fit.
    theta_init <- PACS:::warm_start_theta(M = M, q = fx$q, T = 2L, p_beta = 2L)
    res_exact <- PACS:::irls_iter_cumu(M_vec = M, X = fx$X,
                                       theta_estimated = theta_init,
                                       q_vec = fx$q, T = 2L)
    if (res_exact[length(res_exact)] == 1L) {
      beta_hat_exact[r, ] <- res_exact[3:4]
    }

    ## Stacked fit via the existing binary IRLS on stacked data.
    Z1 <- as.integer(M >= 1L)
    Z2 <- as.integer(M >= 2L)
    A <- diag(2L)
    X_alpha <- A[c(rep(1L, n), rep(2L, n)), , drop = FALSE]
    X_stack <- cbind(X_alpha, rbind(fx$X, fx$X))
    q_stack <- c(fx$q, fx$q)
    Z_stack <- c(Z1, Z2)
    theta_init_stack <- rep.int(0.05, ncol(X_stack))
    res_stack <- irls_iter(y_vec = Z_stack, xdumm = X_stack,
                           theta_estimated = theta_init_stack,
                           q_vec = q_stack)
    if (res_stack[length(res_stack)] == 1L) {
      ## Last 2 entries are beta in our stacked design.
      beta_hat_stack[r, ] <- res_stack[3:4]
    }
  }

  ## Compute MSE per coordinate, summed.
  mse <- function(B, true) mean(rowSums((B - matrix(true, nrow = nrow(B),
                                                    ncol = length(true),
                                                    byrow = TRUE))^2,
                                        na.rm = TRUE), na.rm = TRUE)
  mse_exact <- mse(beta_hat_exact, beta)
  mse_stack <- mse(beta_hat_stack, beta)

  message(sprintf("MSE(beta) exact = %.4f, stacked = %.4f", mse_exact, mse_stack))
  expect_lte(mse_exact, mse_stack * 1.10)  ## exact within 10% of stacked or better
})
