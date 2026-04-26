## Under the null (true beta_test = 0), the cumulative-logit LRT p-values
## from compare_models_cumu should be approximately uniform.

test_that("LRT p-values are roughly uniform under the null", {
  skip_on_cran()
  n_rep <- 100L
  n <- 300L
  alpha <- c(0.7, -0.3)
  beta_full <- c(0.4, 0.0)  ## second covariate is the test target; set to 0
  pvals <- numeric(n_rep)
  ok <- logical(n_rep)

  for (r in seq_len(n_rep)) {
    fx <- make_cumu_fixture(n = n, p = 2L, T = 2L,
                            alpha = alpha, beta = beta_full,
                            seed = 1000L + r)
    M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                   q = fx$q, capture = "B", seed = 1500L + r)

    ## Full fit: both betas free.
    th_init <- PACS:::warm_start_theta(M = M, q = fx$q, T = 2L, p_beta = 2L)
    res_full <- PACS:::irls_iter_cumu(M_vec = M, X = fx$X,
                                      theta_estimated = th_init,
                                      q_vec = fx$q, T = 2L)
    ## Null fit: beta_2 (last entry) held at 0.
    th_init_null <- th_init
    th_init_null[4L] <- 0
    res_null <- PACS:::irls_iter_cumu_null(M_vec = M, X = fx$X,
                                           theta_estimated = th_init_null,
                                           hold_zero = 4L,
                                           q_vec = fx$q, T = 2L)

    if (res_full[length(res_full)] != 1L || res_null[length(res_null)] != 1L) {
      next
    }

    th_full <- res_full[seq_len(length(res_full) - 1L)]
    th_null <- res_null[seq_len(length(res_null) - 1L)]

    pv <- PACS:::compare_models_cumu(
      x_full = fx$X,
      theta_estimated_full = matrix(th_full, ncol = 1L),
      theta_estimated_null = matrix(th_null, ncol = 1L),
      q_vec = fx$q,
      c_by_r = matrix(M, ncol = 1L),
      T = 2L, df_test = 1L, mc.cores = 1L
    )
    if (is.finite(pv) && pv >= 0 && pv <= 1) {
      pvals[r] <- pv
      ok[r] <- TRUE
    }
  }

  pv_ok <- pvals[ok]
  expect_gt(length(pv_ok), 0.7 * n_rep)  ## most replicates produced a p-value

  ## Kolmogorov-Smirnov test for uniformity. Allow some slack: if KS p > 0.01
  ## we consider the calibration acceptable for this small simulation.
  ks_p <- suppressWarnings(stats::ks.test(pv_ok, "punif")$p.value)
  message(sprintf("LRT calibration: n=%d valid p-values, KS p=%.3f",
                  length(pv_ok), ks_p))
  expect_gt(ks_p, 0.01)
})
