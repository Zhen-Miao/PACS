## At q = 1 (no capture-rate layer), the cumulative-logit model reduces to
## standard proportional-odds logistic regression. Validate against
## ordinal::clm() which gives the unpenalized MLE.
##
## Two regimes:
##   - Large n (n=2000): Firth correction is negligible, so our Firth-penalized
##     estimates should nearly match clm().
##   - Small n (n=80): Firth correction is non-negligible, so estimates differ,
##     but our *unpenalized* log-likelihood evaluated at the clm() MLE should
##     match clm()'s logLik exactly.

test_that("at q=1 large n, Firth estimates match clm() closely", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  set.seed(42L)
  n <- 2000L
  T <- 2L
  p <- 2L
  alpha <- c(0.8, -0.4)
  beta <- c(0.6, -0.3)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  q <- rep(1.0, n)

  M <- PACS:::simulate_cumu_pacs(X = X, alpha = alpha, beta = beta,
                                 q = q, capture = "B", seed = 100L)

  ## --- Our exact path (Firth-penalized IRLS) ---
  theta_init <- PACS:::warm_start_theta(M = M, q = q, T = T, p_beta = p)
  res <- PACS:::irls_iter_cumu(M_vec = M, X = X, theta_estimated = theta_init,
                               q_vec = q, T = T)
  expect_equal(res[length(res)], 1)  ## converged
  th <- res[seq_len(length(res) - 1L)]
  alpha_hat <- PACS:::alpha_from_atilde(th[1:T])
  beta_hat <- th[(T + 1):(T + p)]

  ## --- ordinal::clm() ---
  Y_factor <- factor(M, levels = 0:T, ordered = TRUE)
  df <- data.frame(Y = Y_factor, x1 = X[, 1], x2 = X[, 2])
  fit_clm <- ordinal::clm(Y ~ x1 + x2, data = df, link = "logit")

  ## clm() parameterises thresholds as α^clm_t = -α_t (opposite sign convention:
  ## clm uses Pr(Y <= k) = σ(θ_k - x'β), we use Pr(Y >= t) = σ(α_t + x'β)).
  alpha_clm <- -as.numeric(fit_clm$alpha)
  beta_clm <- as.numeric(fit_clm$beta)

  ## At n=2000, Firth shrinkage is O(1/n) ≈ 5e-4, so 0.05 tolerance is generous.
  expect_equal(alpha_hat, alpha_clm, tolerance = 0.05,
               label = "alpha vs clm at n=2000")
  expect_equal(beta_hat, beta_clm, tolerance = 0.05,
               label = "beta vs clm at n=2000")

  message(sprintf("n=2000: alpha PACS=[%.4f, %.4f] clm=[%.4f, %.4f]",
                  alpha_hat[1], alpha_hat[2], alpha_clm[1], alpha_clm[2]))
  message(sprintf("n=2000: beta  PACS=[%.4f, %.4f] clm=[%.4f, %.4f]",
                  beta_hat[1], beta_hat[2], beta_clm[1], beta_clm[2]))
})


test_that("at q=1 small n, unpenalized logLik matches clm() exactly", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  set.seed(77L)
  n <- 80L
  T <- 2L
  p <- 2L
  alpha <- c(0.8, -0.4)
  beta <- c(0.6, -0.3)
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  q <- rep(1.0, n)

  M <- PACS:::simulate_cumu_pacs(X = X, alpha = alpha, beta = beta,
                                 q = q, capture = "B", seed = 200L)

  ## --- ordinal::clm() ---
  Y_factor <- factor(M, levels = 0:T, ordered = TRUE)
  df <- data.frame(Y = Y_factor, x1 = X[, 1], x2 = X[, 2])
  fit_clm <- ordinal::clm(Y ~ x1 + x2, data = df, link = "logit")

  ## Map clm estimates to our parameterisation.
  alpha_clm <- -as.numeric(fit_clm$alpha)
  beta_clm <- as.numeric(fit_clm$beta)

  ## Build theta in our (atilde, beta) space.
  atilde_clm <- numeric(T)
  atilde_clm[1] <- alpha_clm[1]
  for (t in 2:T) {
    atilde_clm[t] <- log(alpha_clm[t - 1] - alpha_clm[t])
  }
  theta_clm <- c(atilde_clm, beta_clm)

  ## Our unpenalized log-likelihood at the clm MLE should match clm's logLik.
  ## At q=1 our loss_fun_cumu drops the log(q_i) term for m>=1, but q=1 =>
  ## log(1)=0, so no correction needed.
  ll_ours <- PACS:::loss_fun_cumu(theta = theta_clm, X = X, M = M, q = q, T = T)
  ll_clm <- as.numeric(logLik(fit_clm))

  expect_equal(ll_ours, ll_clm, tolerance = 1e-6,
               label = "unpenalized logLik at clm MLE, small n")

  message(sprintf("n=80: logLik ours=%.6f, clm=%.6f, diff=%.2e",
                  ll_ours, ll_clm, abs(ll_ours - ll_clm)))
})


test_that("at q=1 small n, score at clm MLE is zero", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  set.seed(88L)
  n <- 80L
  T <- 2L
  p <- 1L
  alpha <- c(0.5, -0.5)
  beta <- 0.4
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- "x1"
  q <- rep(1.0, n)

  M <- PACS:::simulate_cumu_pacs(X = X, alpha = alpha, beta = beta,
                                 q = q, capture = "B", seed = 300L)

  Y_factor <- factor(M, levels = 0:T, ordered = TRUE)
  df <- data.frame(Y = Y_factor, x1 = X[, 1])
  fit_clm <- ordinal::clm(Y ~ x1, data = df, link = "logit")

  alpha_clm <- -as.numeric(fit_clm$alpha)
  beta_clm <- as.numeric(fit_clm$beta)

  atilde_clm <- numeric(T)
  atilde_clm[1] <- alpha_clm[1]
  for (t in 2:T) {
    atilde_clm[t] <- log(alpha_clm[t - 1] - alpha_clm[t])
  }
  theta_clm <- c(atilde_clm, beta_clm)

  grad <- PACS:::loss_gradient_cumu(theta = theta_clm, X = X, M = M, q = q, T = T)
  expect_equal(as.numeric(grad), rep(0, length(grad)), tolerance = 1e-4,
               label = "score at clm MLE should be ~0")

  message(sprintf("n=80: max |score| at clm MLE = %.2e", max(abs(grad))))
})


test_that("at q=1 Firth reduces bias vs clm MLE at small n (simulation)", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  n_rep <- 40L
  n <- 80L
  T <- 2L
  p <- 1L
  alpha <- c(0.8, -0.4)
  beta <- 0.5

  beta_firth <- numeric(n_rep)
  beta_clm <- numeric(n_rep)
  ok_firth <- logical(n_rep)
  ok_clm <- logical(n_rep)

  for (r in seq_len(n_rep)) {
    set.seed(400L + r)
    X <- matrix(rnorm(n * p), nrow = n, ncol = p)
    colnames(X) <- "x1"
    q <- rep(1.0, n)
    M <- PACS:::simulate_cumu_pacs(X = X, alpha = alpha, beta = beta,
                                   q = q, capture = "B", seed = 500L + r)

    ## Firth (our code)
    theta_init <- PACS:::warm_start_theta(M = M, q = q, T = T, p_beta = p)
    res <- PACS:::irls_iter_cumu(M_vec = M, X = X, theta_estimated = theta_init,
                                 q_vec = q, T = T)
    if (res[length(res)] == 1L) {
      beta_firth[r] <- res[T + 1L]
      ok_firth[r] <- TRUE
    }

    ## clm (unpenalized MLE)
    Y_factor <- factor(M, levels = 0:T, ordered = TRUE)
    df <- data.frame(Y = Y_factor, x1 = X[, 1])
    fit <- try(ordinal::clm(Y ~ x1, data = df, link = "logit"), silent = TRUE)
    if (!inherits(fit, "try-error") && fit$convergence$code == 0L) {
      beta_clm[r] <- as.numeric(fit$beta)
      ok_clm[r] <- TRUE
    }
  }

  both_ok <- ok_firth & ok_clm
  expect_gt(sum(both_ok), 0.7 * n_rep)

  bias_firth <- mean(beta_firth[both_ok]) - beta
  bias_clm <- mean(beta_clm[both_ok]) - beta

  message(sprintf("n=80, %d reps: bias(Firth)=%.4f, bias(clm MLE)=%.4f",
                  sum(both_ok), bias_firth, bias_clm))

  ## Firth should have smaller absolute bias than unpenalized MLE.
  expect_lte(abs(bias_firth), abs(bias_clm) + 0.02)
})
