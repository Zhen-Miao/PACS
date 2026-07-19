## End-to-end smoke test for pacs_test_cumu(method = "exact").
## Verifies the public API wiring: formula parsing, intercept stripping,
## parameter-of-interest identification, warm start, fitting, and LRT
## all work together correctly.

test_that("pacs_test_cumu exact returns correct structure", {
  set.seed(99L)
  n <- 200L
  n_peaks <- 5L
  T <- 2L
  alpha <- c(0.6, -0.3)
  beta <- 0.5

  group <- factor(c(rep("A", n %/% 2), rep("B", n - n %/% 2)))
  X_raw <- model.matrix(~ group)[, -1, drop = FALSE]
  q <- runif(n, 0.4, 0.9)

  ## Simulate a small pic_matrix (peaks x cells) with count values in 0..T.
  pic <- matrix(0L, nrow = n_peaks, ncol = n)
  rownames(pic) <- paste0("peak", seq_len(n_peaks))
  for (j in seq_len(n_peaks)) {
    pic[j, ] <- PACS:::simulate_cumu_pacs(
      X = X_raw, alpha = alpha, beta = beta,
      q = q, capture = "B", seed = 900L + j
    )
  }

  meta <- data.frame(group = group)

  result <- pacs_test_cumu(
    covariate_meta.data = meta,
    formula_full = ~ group,
    formula_null = ~ 1,
    pic_matrix = pic,
    max_T = T,
    cap_rates = q,
    n_cores = 1L,
    method = "exact"
  )

  ## --- Structure checks ---
  expect_type(result, "list")
  expect_named(result, c("pacs_converged", "pacs_p_val"))
  expect_length(result$pacs_p_val, n_peaks)
  expect_true(all(names(result$pacs_p_val) == rownames(pic)))

  ## p-values should be in [0, 1]
  expect_true(all(result$pacs_p_val >= 0 & result$pacs_p_val <= 1,
                  na.rm = TRUE))

  ## Convergence vector: 2 * n_peaks (null + full)
  expect_length(result$pacs_converged, 2L * n_peaks)

  message(sprintf("API smoke: %d peaks, p-values = [%s]",
                  n_peaks,
                  paste(sprintf("%.4f", result$pacs_p_val), collapse = ", ")))
})


test_that("pacs_test_cumu exact matches internal functions", {
  skip_on_cran()
  set.seed(55L)
  n <- 150L
  n_peaks <- 3L
  T <- 2L
  alpha <- c(0.7, -0.4)
  beta <- 0.3

  group <- factor(c(rep("ctrl", n %/% 2), rep("treat", n - n %/% 2)))
  X_raw <- model.matrix(~ group)[, -1, drop = FALSE]
  colnames(X_raw) <- "grouptreat"
  q <- runif(n, 0.5, 0.8)

  pic <- matrix(0L, nrow = n_peaks, ncol = n)
  rownames(pic) <- paste0("peak", seq_len(n_peaks))
  for (j in seq_len(n_peaks)) {
    pic[j, ] <- PACS:::simulate_cumu_pacs(
      X = X_raw, alpha = alpha, beta = beta,
      q = q, capture = "B", seed = 1100L + j
    )
  }

  meta <- data.frame(group = group)

  ## Public API
  res_api <- pacs_test_cumu(
    covariate_meta.data = meta,
    formula_full = ~ group,
    formula_null = ~ 1,
    pic_matrix = pic,
    max_T = T,
    cap_rates = q,
    n_cores = 1L,
    method = "exact"
  )

  ## Internal path: fit each peak manually
  pv_manual <- numeric(n_peaks)
  for (j in seq_len(n_peaks)) {
    M_j <- pic[j, ]
    th_init <- PACS:::warm_start_theta(M = M_j, q = q, T = T, p_beta = 1L)
    res_full <- PACS:::irls_iter_cumu(
      M_vec = M_j, X = X_raw, theta_estimated = th_init, q_vec = q, T = T
    )
    th_init_null <- th_init
    th_init_null[3L] <- 0
    res_null <- PACS:::irls_iter_cumu_null(
      M_vec = M_j, X = X_raw, theta_estimated = th_init_null,
      hold_zero = 3L, q_vec = q, T = T
    )
    if (res_full[length(res_full)] == 1L &&
        res_null[length(res_null)] == 1L) {
      th_f <- res_full[1:3]
      th_n <- res_null[1:3]
      pv_manual[j] <- PACS:::compare_models_cumu(
        x_full = X_raw,
        theta_estimated_full = matrix(th_f, ncol = 1L),
        theta_estimated_null = matrix(th_n, ncol = 1L),
        q_vec = q, c_by_r = matrix(M_j, ncol = 1L),
        T = T, df_test = 1L, mc.cores = 1L
      )
    } else {
      pv_manual[j] <- NA_real_
    }
  }

  ## The public API pools all peaks for a shared warm start, so estimates
  ## can differ slightly. But p-values should be close when both converge.
  both_finite <- is.finite(pv_manual) & is.finite(res_api$pacs_p_val)
  expect_gt(sum(both_finite), 0)

  ## Log the comparison rather than assert exact equality (warm-start
  ## differences can propagate), but check they're in the same ballpark.
  for (j in which(both_finite)) {
    message(sprintf("peak%d: API p=%.4f, manual p=%.4f",
                    j, res_api$pacs_p_val[j], pv_manual[j]))
  }

  ## Both should agree on which peaks are significant at a liberal threshold.
  sig_api <- unname(res_api$pacs_p_val[both_finite] < 0.10)
  sig_manual <- unname(pv_manual[both_finite] < 0.10)
  expect_equal(sig_api, sig_manual,
               label = "API and manual agree on significance calls")
})


test_that("pacs_test_cumu exact handles multiple covariates", {
  set.seed(33L)
  n <- 200L
  n_peaks <- 3L
  T <- 2L

  group <- factor(c(rep("A", n %/% 2), rep("B", n - n %/% 2)))
  batch <- rnorm(n)
  meta <- data.frame(group = group, batch = batch)

  X_full <- model.matrix(~ group + batch, data = meta)
  X_raw <- X_full[, -1, drop = FALSE]
  q <- runif(n, 0.5, 0.8)

  ## Simulate with group effect, no batch effect
  alpha <- c(0.6, -0.3)
  beta <- c(0.4, 0.0)  ## groupB effect, no batch effect
  pic <- matrix(0L, nrow = n_peaks, ncol = n)
  rownames(pic) <- paste0("peak", seq_len(n_peaks))
  for (j in seq_len(n_peaks)) {
    pic[j, ] <- PACS:::simulate_cumu_pacs(
      X = X_raw, alpha = alpha, beta = beta,
      q = q, capture = "B", seed = 1200L + j
    )
  }

  ## Test group effect adjusting for batch
  result <- pacs_test_cumu(
    covariate_meta.data = meta,
    formula_full = ~ group + batch,
    formula_null = ~ batch,
    pic_matrix = pic,
    max_T = T,
    cap_rates = q,
    n_cores = 1L,
    method = "exact"
  )

  expect_length(result$pacs_p_val, n_peaks)
  expect_true(all(result$pacs_p_val >= 0 & result$pacs_p_val <= 1,
                  na.rm = TRUE))

  message(sprintf("Multi-covariate: p-values = [%s]",
                  paste(sprintf("%.4f", result$pacs_p_val), collapse = ", ")))
})
