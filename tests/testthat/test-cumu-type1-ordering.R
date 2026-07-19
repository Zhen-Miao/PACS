## Under the null (true beta_test = 0), the stacked method inflates type-1
## error because it treats the nested indicators Z_{i1} >= Z_{i2} as
## independent, double-counting each cell. The binary method (using only
## M >= 1) is conservative because it discards threshold information.
##
## Expected ordering of rejection rates at any nominal level alpha:
##   binary <= exact <= stacked
##
## Equivalently, the median p-value should satisfy:
##   stacked <= exact <= binary
##
## We verify this across multiple sample sizes (n = 100, 300, 800).

test_that("type-1 error ordering: binary <= exact <= stacked across sample sizes", {
  skip_on_cran()

  alpha_true <- c(0.7, -0.3)
  beta_nuisance <- 0.4
  beta_test <- 0.0  ## null hypothesis
  T <- 2L
  n_rep <- 150L
  nominal_alpha <- 0.05

  sample_sizes <- c(100L, 300L, 800L)

  for (n in sample_sizes) {
    pv_exact <- numeric(n_rep)
    pv_stack <- numeric(n_rep)
    pv_binary <- numeric(n_rep)
    ok_exact <- logical(n_rep)
    ok_stack <- logical(n_rep)
    ok_binary <- logical(n_rep)

    for (r in seq_len(n_rep)) {
      set.seed(2000L + r + n)
      X <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)
      colnames(X) <- c("x1", "x2")
      q <- runif(n, 0.4, 0.9)

      M <- PACS:::simulate_cumu_pacs(
        X = X, alpha = alpha_true,
        beta = c(beta_nuisance, beta_test),
        q = q, capture = "B", seed = 3000L + r + n
      )

      ## ---- Exact cumulative-logit path ----
      th_init <- PACS:::warm_start_theta(M = M, q = q, T = T, p_beta = 2L)
      res_full <- PACS:::irls_iter_cumu(
        M_vec = M, X = X, theta_estimated = th_init, q_vec = q, T = T
      )
      th_init_null <- th_init
      th_init_null[4L] <- 0
      res_null <- PACS:::irls_iter_cumu_null(
        M_vec = M, X = X, theta_estimated = th_init_null,
        hold_zero = 4L, q_vec = q, T = T
      )
      if (res_full[length(res_full)] == 1L &&
          res_null[length(res_null)] == 1L) {
        th_f <- res_full[1:4]
        th_n <- res_null[1:4]
        pv <- PACS:::compare_models_cumu(
          x_full = X,
          theta_estimated_full = matrix(th_f, ncol = 1L),
          theta_estimated_null = matrix(th_n, ncol = 1L),
          q_vec = q, c_by_r = matrix(M, ncol = 1L),
          T = T, df_test = 1L, mc.cores = 1L
        )
        if (is.finite(pv) && pv >= 0 && pv <= 1) {
          pv_exact[r] <- pv
          ok_exact[r] <- TRUE
        }
      }

      ## ---- Stacked path (binary IRLS on duplicated rows) ----
      ## Uses the same full-design-matrix + hold_zero approach as pacs_test_logit.
      Z1 <- as.integer(M >= 1L)
      Z2 <- as.integer(M >= 2L)
      A <- diag(T)
      X_alpha <- A[c(rep(1L, n), rep(2L, n)), , drop = FALSE]
      X_stack <- cbind(X_alpha, rbind(X, X))
      q_stack <- c(q, q)
      Z_stack <- c(Z1, Z2)

      ## Full stacked fit (all 4 params free: alpha1, alpha2, beta_x1, beta_x2)
      theta_init_s <- rep.int(0.05, ncol(X_stack))
      res_s_full <- irls_iter(y_vec = Z_stack, xdumm = X_stack,
                              theta_estimated = theta_init_s, q_vec = q_stack)
      ## Null stacked fit: hold beta_x2 (col 4) at zero, same design matrix
      theta_init_s_null <- rep.int(0.05, ncol(X_stack))
      theta_init_s_null[4L] <- 0
      res_s_null <- irls_iter_null(y_vec = Z_stack, xdumm = X_stack,
                                   theta_estimated = theta_init_s_null,
                                   hold_zero = 4L, q_vec = q_stack)
      if (res_s_full[length(res_s_full)] == 1L &&
          res_s_null[length(res_s_null)] == 1L) {
        th_sf <- res_s_full[seq_len(length(res_s_full) - 1L)]
        th_sn <- res_s_null[seq_len(length(res_s_null) - 1L)]
        pv_s <- compare_models(
          x_full = X_stack,
          theta_estimated_full = matrix(th_sf, ncol = 1L),
          x_null = X_stack,
          theta_estimated_null = matrix(th_sn, ncol = 1L),
          q_vec = q_stack, c_by_r = matrix(Z_stack, ncol = 1L),
          df_test = 1L, mc.cores = 1L
        )
        if (is.finite(pv_s)) {
          pv_stack[r] <- pv_s
          ok_stack[r] <- TRUE
        }
      }

      ## ---- Binary path (only Z1 = 1[M >= 1]) ----
      ## Uses full design matrix + hold_zero, matching pacs_test_logit convention.
      Z_bin <- as.integer(M >= 1L)
      X_bin <- cbind(intercept = 1, X)  ## 3 cols: intercept, x1, x2

      ## Full binary fit (intercept + x1 + x2)
      theta_init_b <- rep.int(0.05, ncol(X_bin))
      res_b_full <- irls_iter(y_vec = Z_bin, xdumm = X_bin,
                              theta_estimated = theta_init_b, q_vec = q)
      ## Null binary fit: hold x2 (col 3) at zero, same design matrix
      theta_init_b_null <- rep.int(0.05, ncol(X_bin))
      theta_init_b_null[3L] <- 0
      res_b_null <- irls_iter_null(y_vec = Z_bin, xdumm = X_bin,
                                   theta_estimated = theta_init_b_null,
                                   hold_zero = 3L, q_vec = q)
      if (res_b_full[length(res_b_full)] == 1L &&
          res_b_null[length(res_b_null)] == 1L) {
        th_bf <- res_b_full[seq_len(length(res_b_full) - 1L)]
        th_bn <- res_b_null[seq_len(length(res_b_null) - 1L)]
        pv_b <- compare_models(
          x_full = X_bin,
          theta_estimated_full = matrix(th_bf, ncol = 1L),
          x_null = X_bin,
          theta_estimated_null = matrix(th_bn, ncol = 1L),
          q_vec = q, c_by_r = matrix(Z_bin, ncol = 1L),
          df_test = 1L, mc.cores = 1L
        )
        if (is.finite(pv_b)) {
          pv_binary[r] <- pv_b
          ok_binary[r] <- TRUE
        }
      }
    }

    all_ok <- ok_exact & ok_stack & ok_binary
    n_ok <- sum(all_ok)
    expect_gt(n_ok, 0.6 * n_rep,
              label = sprintf("convergence rate at n=%d", n))

    rej_exact <- mean(pv_exact[all_ok] < nominal_alpha)
    rej_stack <- mean(pv_stack[all_ok] < nominal_alpha)
    rej_binary <- mean(pv_binary[all_ok] < nominal_alpha)

    med_exact <- median(pv_exact[all_ok])
    med_stack <- median(pv_stack[all_ok])
    med_binary <- median(pv_binary[all_ok])

    message(sprintf(
      paste0("n=%d (%d reps): rejection rate at alpha=%.2f: ",
             "binary=%.3f, exact=%.3f, stacked=%.3f  |  ",
             "median p: binary=%.3f, exact=%.3f, stacked=%.3f"),
      n, n_ok, nominal_alpha,
      rej_binary, rej_exact, rej_stack,
      med_binary, med_exact, med_stack
    ))

    ## Core assertion: stacked rejection rate >= exact (with slack for noise).
    expect_gte(rej_stack + 0.03, rej_exact,
               label = sprintf("stacked >= exact rejection at n=%d", n))
    ## Exact rejection rate >= binary (with slack).
    expect_gte(rej_exact + 0.05, rej_binary,
               label = sprintf("exact >= binary rejection at n=%d", n))

    ## Median p-value ordering: binary >= exact >= stacked (with slack).
    expect_gte(med_binary + 0.10, med_exact,
               label = sprintf("binary median p >= exact at n=%d", n))
    expect_gte(med_exact + 0.10, med_stack,
               label = sprintf("exact median p >= stacked at n=%d", n))
  }
})
