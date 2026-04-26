## When two groups have different capture rates but identical accessibility
## (beta_group = 0), clm() — which has no capture-rate model — confounds
## "captured more" with "more accessible" and produces a spurious group
## effect. PACS corrects for this via the q_i layer.
##
## Test design: vary (q_A, q_B) from equal to highly divergent. At each
## combination, simulate n_rep replicates under the null and measure:
##   1. Mean |beta_group| estimate (bias)
##   2. Rejection rate at nominal 5%
##
## Expected trends:
##   - PACS: beta_group ≈ 0 and rejection ≈ 5% regardless of q gap
##   - clm:  bias and rejection grow monotonically with |q_A - q_B|
##   - At equal q: PACS ≈ clm (both correct)

test_that("PACS stays calibrated while clm() degrades as capture-rate gap widens", {
  skip_on_cran()
  skip_if_not_installed("ordinal")

  n <- 400L
  T <- 2L
  alpha_true <- c(0.7, -0.3)
  beta_group_true <- 0.0  ## null: no accessibility difference
  n_rep <- 80L
  nominal_alpha <- 0.05

  ## Capture-rate scenarios: (q_groupA, q_groupB), ordered by increasing gap.
  scenarios <- list(
    c(0.70, 0.70),
    c(0.60, 0.80),
    c(0.50, 0.80),
    c(0.40, 0.80),
    c(0.30, 0.80)
  )

  rej_pacs <- numeric(length(scenarios))
  rej_clm <- numeric(length(scenarios))
  bias_pacs <- numeric(length(scenarios))
  bias_clm <- numeric(length(scenarios))

  for (s in seq_along(scenarios)) {
    q_A <- scenarios[[s]][1]
    q_B <- scenarios[[s]][2]

    pv_pacs <- numeric(n_rep)
    pv_clm <- numeric(n_rep)
    beta_hat_pacs <- numeric(n_rep)
    beta_hat_clm <- numeric(n_rep)
    ok_pacs <- logical(n_rep)
    ok_clm <- logical(n_rep)

    for (r in seq_len(n_rep)) {
      set.seed(5000L + r + s * 1000L)
      n_A <- n %/% 2L
      n_B <- n - n_A

      ## Group indicator (the covariate of interest)
      group <- c(rep(0L, n_A), rep(1L, n_B))
      X <- matrix(group, ncol = 1L)
      colnames(X) <- "group"

      ## Per-cell capture rates
      q <- c(rep(q_A, n_A), rep(q_B, n_B))

      ## Simulate under null: same alpha for both groups, beta_group = 0
      M <- PACS:::simulate_cumu_pacs(
        X = X, alpha = alpha_true,
        beta = beta_group_true,
        q = q, capture = "B", seed = 6000L + r + s * 1000L
      )

      ## ---- PACS exact cumulative-logit (knows about q) ----
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
        beta_hat_pacs[r] <- th_f[3L]
        pv <- PACS:::compare_models_cumu(
          x_full = X,
          theta_estimated_full = matrix(th_f, ncol = 1L),
          theta_estimated_null = matrix(th_n, ncol = 1L),
          q_vec = q, c_by_r = matrix(M, ncol = 1L),
          T = T, df_test = 1L, mc.cores = 1L
        )
        if (is.finite(pv) && pv >= 0 && pv <= 1) {
          pv_pacs[r] <- pv
          ok_pacs[r] <- TRUE
        }
      }

      ## ---- ordinal::clm() (ignores capture rate) ----
      Y_factor <- factor(M, levels = 0:T, ordered = TRUE)
      df <- data.frame(Y = Y_factor, group = group)
      fit <- try(ordinal::clm(Y ~ group, data = df, link = "logit"),
                 silent = TRUE)
      if (!inherits(fit, "try-error") && fit$convergence$code == 0L) {
        beta_hat_clm[r] <- as.numeric(fit$beta)
        ## Wald test p-value from clm summary
        s_fit <- summary(fit)
        pv_clm[r] <- s_fit$coefficients["group", "Pr(>|z|)"]
        ok_clm[r] <- TRUE
      }
    }

    both_ok <- ok_pacs & ok_clm
    n_ok <- sum(both_ok)

    rej_pacs[s] <- mean(pv_pacs[both_ok] < nominal_alpha)
    rej_clm[s] <- mean(pv_clm[both_ok] < nominal_alpha)
    bias_pacs[s] <- mean(abs(beta_hat_pacs[both_ok]))
    bias_clm[s] <- mean(abs(beta_hat_clm[both_ok]))

    message(sprintf(
      "q=(%0.2f,%0.2f) gap=%.2f [%d reps]: rej PACS=%.3f clm=%.3f | mean|beta| PACS=%.4f clm=%.4f",
      q_A, q_B, abs(q_A - q_B), n_ok,
      rej_pacs[s], rej_clm[s], bias_pacs[s], bias_clm[s]
    ))
  }

  ## ---- Assertions ----

  ## 1. At equal capture rates (scenario 1), PACS and clm should be similar:
  ##    both well-calibrated, both low bias.
  expect_lt(rej_pacs[1], 0.15, label = "PACS rejection at equal q")
  expect_lt(rej_clm[1], 0.15, label = "clm rejection at equal q")
  expect_lt(abs(bias_pacs[1] - bias_clm[1]), 0.10,
            label = "bias gap at equal q")

  ## 2. PACS stays calibrated across all scenarios.
  expect_lt(max(rej_pacs), 0.15,
            label = "PACS rejection stays controlled")

  ## 3. clm rejection rate at the widest gap should exceed PACS.
  expect_gt(rej_clm[length(scenarios)], rej_pacs[length(scenarios)],
            label = "clm rejects more than PACS at wide q gap")

  ## 4. clm bias should increase monotonically with the gap (allow ties).
  ##    Check that the widest-gap bias exceeds the equal-q bias.
  expect_gt(bias_clm[length(scenarios)], bias_clm[1] + 0.01,
            label = "clm bias grows with q gap")

  ## 5. PACS bias should stay roughly flat (contrast: clm bias spans >1.0).
  expect_lt(max(bias_pacs) - min(bias_pacs), 0.15,
            label = "PACS bias stable across q gaps")

  ## Print summary table for review.
  message("\n--- Summary ---")
  message("q_gap  | rej_PACS  rej_clm  | bias_PACS  bias_clm")
  for (s in seq_along(scenarios)) {
    message(sprintf("  %.2f |   %.3f    %.3f   |   %.4f    %.4f",
                    abs(scenarios[[s]][1] - scenarios[[s]][2]),
                    rej_pacs[s], rej_clm[s], bias_pacs[s], bias_clm[s]))
  }
})
