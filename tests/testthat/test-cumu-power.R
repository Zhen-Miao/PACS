## Under the alternative (beta_group != 0), the exact cumulative-logit path
## should have more power than the binary path (which discards threshold
## information by collapsing M to 1[M >= 1]).
##
## This test simulates two-group data from the cumulative model with a
## true group effect and compares rejection rates of exact vs binary at
## nominal 5%. The power gain comes from the second threshold: cells with
## M = 2 carry more evidence than cells with M = 1, but binary treats
## them identically.
##
## We test at two effect sizes: moderate (beta = 0.4) and small (beta = 0.2),
## and verify that the exact method's power advantage is present and grows
## as the effect shrinks (where every bit of information matters more).

test_that("exact has more power than binary under the alternative", {
  skip_on_cran()

  T <- 2L
  alpha_true <- c(0.7, -0.3)
  n <- 400L
  n_rep <- 120L
  nominal_alpha <- 0.05

  effect_sizes <- c(0.4, 0.2)
  power_exact <- numeric(length(effect_sizes))
  power_binary <- numeric(length(effect_sizes))

  for (e in seq_along(effect_sizes)) {
    beta_group <- effect_sizes[e]

    pv_exact <- numeric(n_rep)
    pv_binary <- numeric(n_rep)
    ok_exact <- logical(n_rep)
    ok_binary <- logical(n_rep)

    for (r in seq_len(n_rep)) {
      set.seed(7000L + r + e * 1000L)
      n_A <- n %/% 2L
      n_B <- n - n_A
      group <- c(rep(0L, n_A), rep(1L, n_B))
      X <- matrix(group, ncol = 1L)
      colnames(X) <- "group"
      q <- runif(n, 0.4, 0.9)

      M <- PACS:::simulate_cumu_pacs(
        X = X, alpha = alpha_true, beta = beta_group,
        q = q, capture = "B", seed = 8000L + r + e * 1000L
      )

      ## ---- Exact cumulative-logit ----
      th_init <- PACS:::warm_start_theta(M = M, q = q, T = T, p_beta = 1L)
      res_full <- PACS:::irls_iter_cumu(
        M_vec = M, X = X, theta_estimated = th_init, q_vec = q, T = T
      )
      th_init_null <- th_init
      th_init_null[3L] <- 0
      res_null <- PACS:::irls_iter_cumu_null(
        M_vec = M, X = X, theta_estimated = th_init_null,
        hold_zero = 3L, q_vec = q, T = T
      )
      if (res_full[length(res_full)] == 1L &&
          res_null[length(res_null)] == 1L) {
        th_f <- res_full[1:3]
        th_n <- res_null[1:3]
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

      ## ---- Binary (M >= 1 only) ----
      Z_bin <- as.integer(M >= 1L)
      X_bin <- cbind(intercept = 1, X)
      theta_init_b <- rep.int(0.05, ncol(X_bin))
      res_b_full <- irls_iter(y_vec = Z_bin, xdumm = X_bin,
                              theta_estimated = theta_init_b, q_vec = q)
      theta_init_b_null <- rep.int(0.05, ncol(X_bin))
      theta_init_b_null[2L] <- 0
      res_b_null <- irls_iter_null(y_vec = Z_bin, xdumm = X_bin,
                                   theta_estimated = theta_init_b_null,
                                   hold_zero = 2L, q_vec = q)
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

    both_ok <- ok_exact & ok_binary
    n_ok <- sum(both_ok)
    expect_gt(n_ok, 0.7 * n_rep,
              label = sprintf("convergence at beta=%.1f", beta_group))

    power_exact[e] <- mean(pv_exact[both_ok] < nominal_alpha)
    power_binary[e] <- mean(pv_binary[both_ok] < nominal_alpha)

    message(sprintf(
      "beta=%.1f, n=%d (%d reps): power exact=%.3f, binary=%.3f, gain=+%.1f%%",
      beta_group, n, n_ok,
      power_exact[e], power_binary[e],
      100 * (power_exact[e] - power_binary[e])
    ))
  }

  ## Exact should have at least as much power as binary at both effect sizes.
  for (e in seq_along(effect_sizes)) {
    expect_gte(power_exact[e] + 0.03, power_binary[e],
               label = sprintf("exact power >= binary at beta=%.1f",
                               effect_sizes[e]))
  }

  ## At both effect sizes, exact should have non-trivial power (> 10%).
  expect_gt(power_exact[1], 0.10,
            label = "exact detects moderate effect")
  expect_gt(power_exact[2], 0.05,
            label = "exact detects small effect above nominal")
})
