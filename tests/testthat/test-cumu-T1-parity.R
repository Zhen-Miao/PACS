## At T = 1, the cumulative-logit functions in R/param_estimate_cumu.R must
## reduce term-for-term to the binary functions in
## R/param-estimate_logit_get_p_by_t_June.R. Math note §8.

test_that("loss_fun_cumu agrees with loss_fun at T=1", {
  fx <- make_cumu_fixture(n = 50L, p = 2L, T = 1L,
                          alpha = 0.7, beta = c(0.5, -0.2), seed = 11L)
  M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                 q = fx$q, capture = "B", seed = 12L)

  ## Build binary theta (intercept = alpha_1) and X_bin = cbind(1, X).
  theta_cumu <- c(fx$alpha, fx$beta)
  theta_bin <- theta_cumu
  X_bin <- cbind(intercept = 1, fx$X)
  p_bg <- 1 / (1 + exp(-as.numeric(X_bin %*% theta_bin)))

  ll_cumu <- PACS:::loss_fun_cumu(theta = theta_cumu, X = fx$X, M = M,
                                  q = fx$q, T = 1L)
  ll_bin <- loss_fun(p_bg = p_bg, q_vec = fx$q, y_vec = M)

  ## Cumu drops the constant log q_i for m >= 1; binary keeps it.
  expect_equal(ll_cumu, ll_bin - sum(log(fx$q[M == 1L])), tolerance = 1e-10)
})


test_that("loss_gradient_cumu agrees with loss_gradient at T=1", {
  fx <- make_cumu_fixture(n = 50L, p = 2L, T = 1L,
                          alpha = 0.5, beta = c(0.4, -0.3), seed = 21L)
  M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                 q = fx$q, capture = "B", seed = 22L)
  theta_cumu <- c(fx$alpha, fx$beta)
  X_bin <- cbind(intercept = 1, fx$X)
  p_bg <- 1 / (1 + exp(-as.numeric(X_bin %*% theta_cumu)))

  g_cumu <- PACS:::loss_gradient_cumu(theta = theta_cumu, X = fx$X, M = M,
                                      q = fx$q, T = 1L)
  g_bin <- loss_gradient(xdumm = X_bin, p_bg = p_bg, q_vec = fx$q, y_vec = M)
  expect_equal(as.numeric(g_cumu), as.numeric(g_bin), tolerance = 1e-10)
})


test_that("infor_mat_cumu agrees with infor_mat at T=1", {
  fx <- make_cumu_fixture(n = 50L, p = 2L, T = 1L,
                          alpha = 0.3, beta = c(0.2, 0.1), seed = 31L)
  theta_cumu <- c(fx$alpha, fx$beta)
  X_bin <- cbind(intercept = 1, fx$X)
  p_bg <- 1 / (1 + exp(-as.numeric(X_bin %*% theta_cumu)))

  I_cumu <- PACS:::infor_mat_cumu(theta = theta_cumu, X = fx$X, q = fx$q, T = 1L)
  I_bin <- infor_mat(xdumm = X_bin, p_bg = p_bg, q_vec = fx$q)
  expect_equal(as.numeric(I_cumu), as.numeric(I_bin), tolerance = 1e-10)
})


test_that("irls_iter_cumu and irls_iter agree at T=1", {
  fx <- make_cumu_fixture(n = 200L, p = 2L, T = 1L,
                          alpha = 0.4, beta = c(0.5, -0.4), seed = 41L)
  M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                 q = fx$q, capture = "B", seed = 42L)

  theta_init <- rep.int(0.05, 1L + ncol(fx$X))
  X_bin <- cbind(intercept = 1, fx$X)

  res_cumu <- PACS:::irls_iter_cumu(M_vec = M, X = fx$X,
                                    theta_estimated = theta_init,
                                    q_vec = fx$q, T = 1L)
  res_bin <- irls_iter(y_vec = M, xdumm = X_bin,
                       theta_estimated = theta_init,
                       q_vec = fx$q)

  conv_cumu <- res_cumu[length(res_cumu)]
  conv_bin <- res_bin[length(res_bin)]
  expect_true(conv_cumu == 1L,
              info = paste("cumu convergence flag:", conv_cumu))
  expect_true(conv_bin == 1L,
              info = paste("bin convergence flag:", conv_bin))

  th_cumu <- res_cumu[seq_len(length(res_cumu) - 1L)]
  th_bin <- res_bin[seq_len(length(res_bin) - 1L)]
  ## Tolerance loosened from 1e-6 to 1e-4: cumu uses FD on the Firth penalty
  ## derivative (loss_grad_pen_cumu), binary uses an analytic form. Both are
  ## valid optima of the same objective and should agree at FD precision.
  expect_equal(th_cumu, th_bin, tolerance = 1e-4)
})
