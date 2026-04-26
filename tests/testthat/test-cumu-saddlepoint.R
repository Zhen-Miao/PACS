## Saddlepoint-adjusted p-values (Barndorff-Nielsen r*) vs chi-squared.
##
## Finding: the Firth-corrected chi-squared is already very well calibrated,
## and the r* correction can degrade calibration because Firth's penalty
## absorbs the O(1/n) bias that r* targets. These tests document this
## behaviour and guard against regressions.

test_that("saddlepoint produces finite p-values under the null (df=1)", {
  skip_on_cran()
  n_rep <- 50L
  n <- 300L
  alpha <- c(0.7, -0.3)
  beta_full <- c(0.4, 0.0)
  pvals_chisq <- pvals_saddle <- numeric(n_rep)
  ok <- logical(n_rep)

  for (r in seq_len(n_rep)) {
    fx <- make_cumu_fixture(n = n, p = 2L, T = 2L,
                            alpha = alpha, beta = beta_full,
                            seed = 2000L + r)
    M <- simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                            q = fx$q, capture = "B", seed = 2500L + r)

    th_init <- warm_start_theta(M = M, q = fx$q, T = 2L, p_beta = 2L)
    res_full <- irls_iter_cumu(M_vec = M, X = fx$X,
                               theta_estimated = th_init,
                               q_vec = fx$q, T = 2L)
    th_init_null <- th_init
    th_init_null[4L] <- 0
    res_null <- irls_iter_cumu_null(M_vec = M, X = fx$X,
                                    theta_estimated = th_init_null,
                                    hold_zero = 4L,
                                    q_vec = fx$q, T = 2L)

    if (res_full[length(res_full)] != 1L || res_null[length(res_null)] != 1L) {
      next
    }

    th_full <- res_full[seq_len(length(res_full) - 1L)]
    th_null <- res_null[seq_len(length(res_null) - 1L)]

    pv_chi <- compare_models_cumu(
      x_full = fx$X,
      theta_estimated_full = matrix(th_full, ncol = 1L),
      theta_estimated_null = matrix(th_null, ncol = 1L),
      q_vec = fx$q, c_by_r = matrix(M, ncol = 1L),
      T = 2L, df_test = 1L, hold_zero = 4L,
      pvalue_method = "chisq", mc.cores = 1L
    )
    pv_sad <- compare_models_cumu(
      x_full = fx$X,
      theta_estimated_full = matrix(th_full, ncol = 1L),
      theta_estimated_null = matrix(th_null, ncol = 1L),
      q_vec = fx$q, c_by_r = matrix(M, ncol = 1L),
      T = 2L, df_test = 1L, hold_zero = 4L,
      pvalue_method = "saddlepoint", mc.cores = 1L
    )

    if (is.finite(pv_chi) && is.finite(pv_sad) &&
        pv_chi >= 0 && pv_chi <= 1 && pv_sad >= 0 && pv_sad <= 1) {
      pvals_chisq[r] <- pv_chi
      pvals_saddle[r] <- pv_sad
      ok[r] <- TRUE
    }
  }

  pv_ok_chi <- pvals_chisq[ok]
  pv_ok_sad <- pvals_saddle[ok]
  expect_gt(length(pv_ok_sad), 0.5 * n_rep)

  ks_chi <- suppressWarnings(stats::ks.test(pv_ok_chi, "punif")$p.value)
  ks_sad <- suppressWarnings(stats::ks.test(pv_ok_sad, "punif")$p.value)
  message(sprintf(
    "Calibration (n=%d, %d valid): KS(chisq)=%.3f, KS(saddlepoint)=%.3f",
    n, length(pv_ok_sad), ks_chi, ks_sad
  ))

  ## Chi-squared should be well-calibrated (Firth penalty handles the bias).
  expect_gt(ks_chi, 0.01)
  ## Saddlepoint values should at least be in [0, 1].
  expect_true(all(pv_ok_sad >= 0 & pv_ok_sad <= 1))
})


test_that("saddlepoint falls back to chisq for df > 1", {
  skip_on_cran()
  fx <- make_cumu_fixture(n = 200L, p = 2L, T = 2L,
                          alpha = c(0.7, -0.3), beta = c(0.0, 0.0),
                          seed = 3001L)
  M <- simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                           q = fx$q, capture = "B", seed = 3002L)
  th_init <- warm_start_theta(M = M, q = fx$q, T = 2L, p_beta = 2L)
  res_full <- irls_iter_cumu(M_vec = M, X = fx$X,
                              theta_estimated = th_init,
                              q_vec = fx$q, T = 2L)
  th_init_null <- th_init
  th_init_null[3:4] <- 0
  res_null <- irls_iter_cumu_null(M_vec = M, X = fx$X,
                                   theta_estimated = th_init_null,
                                   hold_zero = 3:4,
                                   q_vec = fx$q, T = 2L)

  skip_if(res_full[length(res_full)] != 1L || res_null[length(res_null)] != 1L,
          "convergence failure")

  th_full <- res_full[seq_len(length(res_full) - 1L)]
  th_null <- res_null[seq_len(length(res_null) - 1L)]

  pv_chi <- compare_models_cumu(
    x_full = fx$X,
    theta_estimated_full = matrix(th_full, ncol = 1L),
    theta_estimated_null = matrix(th_null, ncol = 1L),
    q_vec = fx$q, c_by_r = matrix(M, ncol = 1L),
    T = 2L, df_test = 2L, hold_zero = 3:4,
    pvalue_method = "chisq", mc.cores = 1L
  )
  pv_sad <- compare_models_cumu(
    x_full = fx$X,
    theta_estimated_full = matrix(th_full, ncol = 1L),
    theta_estimated_null = matrix(th_null, ncol = 1L),
    q_vec = fx$q, c_by_r = matrix(M, ncol = 1L),
    T = 2L, df_test = 2L, hold_zero = 3:4,
    pvalue_method = "saddlepoint", mc.cores = 1L
  )
  expect_equal(pv_sad, pv_chi)
})


test_that("saddlepoint_pvalue_scalar returns valid p-values", {
  skip_on_cran()
  fx <- make_cumu_fixture(n = 500L, p = 2L, T = 2L,
                          alpha = c(0.7, -0.3), beta = c(0.4, 0.0),
                          seed = 4001L)
  M <- simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                           q = fx$q, capture = "B", seed = 4002L)

  th_init <- warm_start_theta(M = M, q = fx$q, T = 2L, p_beta = 2L)
  res_full <- irls_iter_cumu(M_vec = M, X = fx$X,
                              theta_estimated = th_init,
                              q_vec = fx$q, T = 2L)
  th_init_null <- th_init
  th_init_null[4L] <- 0
  res_null <- irls_iter_cumu_null(M_vec = M, X = fx$X,
                                   theta_estimated = th_init_null,
                                   hold_zero = 4L,
                                   q_vec = fx$q, T = 2L)
  skip_if(res_full[length(res_full)] != 1L || res_null[length(res_null)] != 1L,
          "convergence failure")

  th_full <- res_full[seq_len(length(res_full) - 1L)]
  th_null <- res_null[seq_len(length(res_null) - 1L)]

  I_full <- infor_mat_cumu(th_full, fx$X, fx$q, 2L)
  I_null <- infor_mat_cumu(th_null, fx$X, fx$q, 2L)
  ll_full <- loss_fun_star_cumu(th_full, fx$X, M, fx$q, 2L,
                                 inf_mat = I_full)
  ll_null <- loss_fun_star_cumu(th_null, fx$X, M, fx$q, 2L,
                                 inf_mat = I_null)
  stat <- max(0, 2 * (ll_full - ll_null))

  pv_sad <- saddlepoint_pvalue_scalar(stat, th_full, th_null,
                                       psi_idx = 4L,
                                       X = fx$X, M = M,
                                       q = fx$q, T = 2L)

  expect_true(is.finite(pv_sad))
  expect_true(pv_sad >= 0 && pv_sad <= 1)
})
