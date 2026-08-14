test_that("exact probabilities remain finite near separation", {
  X <- matrix(numeric(), nrow = 1L, ncol = 0L)
  theta <- c(40, log(1e-3))

  link <- PACS:::cumu_link(
    alpha = PACS:::alpha_from_atilde(theta),
    beta = numeric(), X = X
  )
  expect_gt(link$Delta[1L, 1L], 0)

  ll_middle <- PACS:::loss_fun_cumu(
    theta = theta, X = X, M = 1L, q = 1, T = 2L
  )
  ll_zero <- PACS:::loss_fun_cumu(
    theta = theta, X = X, M = 0L, q = 1, T = 2L
  )
  score <- PACS:::loss_gradient_cumu(
    theta = theta, X = X, M = 1L, q = 1, T = 2L
  )
  info <- PACS:::infor_mat_cumu(theta = theta, X = X, q = 1, T = 2L)

  expect_true(is.finite(ll_middle))
  expect_equal(ll_zero, -40, tolerance = 1e-12)
  expect_true(all(is.finite(score)))
  expect_true(all(is.finite(info)))
})


test_that("internal exact helpers reject responses above T clearly", {
  expect_error(
    PACS:::loss_fun_cumu(
      theta = c(0, 0), X = matrix(numeric(), 1L, 0L),
      M = 3L, q = 0.8, T = 2L
    ),
    "M must be in 0:2"
  )
})


test_that("public exact API top-codes counts above max_T", {
  set.seed(121L)
  n <- 240L
  T <- 2L
  group <- factor(rep(c("A", "B"), each = n / 2L))
  meta <- data.frame(group = group)
  X <- model.matrix(~ group, meta)[, -1L, drop = FALSE]
  q <- runif(n, 0.55, 0.9)

  pic_top <- matrix(0L, nrow = 2L, ncol = n)
  rownames(pic_top) <- c("peak1", "peak2")
  for (j in seq_len(nrow(pic_top))) {
    pic_top[j, ] <- PACS:::simulate_cumu_pacs(
      X = X, alpha = c(0.7, -0.2), beta = 0.25,
      q = q, capture = "B", seed = 700L + j
    )
  }
  pic_high <- pic_top
  top_entries <- which(pic_high == T)
  expect_gt(length(top_entries), 0L)
  pic_high[top_entries[seq_len(min(5L, length(top_entries)))]] <- 5L

  fit_top <- suppressWarnings(pacs_test_cumu(
    covariate_meta.data = meta, formula_full = ~ group,
    formula_null = ~ 1, pic_matrix = pic_top, max_T = T,
    cap_rates = q, n_cores = 1L, method = "exact"
  ))
  fit_high <- suppressWarnings(pacs_test_cumu(
    covariate_meta.data = meta, formula_full = ~ group,
    formula_null = ~ 1, pic_matrix = pic_high, max_T = T,
    cap_rates = q, n_cores = 1L, method = "exact"
  ))
  fit_sparse <- suppressWarnings(pacs_test_cumu(
    covariate_meta.data = meta, formula_full = ~ group,
    formula_null = ~ 1,
    pic_matrix = Matrix::Matrix(pic_top, sparse = TRUE), max_T = T,
    cap_rates = q, n_cores = 1L, method = "exact"
  ))

  expect_equal(fit_high, fit_top, tolerance = 1e-12)
  expect_equal(fit_sparse, fit_top, tolerance = 1e-12)
})


test_that("singular exact fits keep singular status", {
  set.seed(91L)
  n <- 80L
  x <- rnorm(n)
  X <- cbind(x, x)
  q <- rep(0.8, n)
  M <- rbinom(n, size = 1L, prob = 0.4)

  fit <- PACS:::irls_iter_cumu(
    M_vec = M, X = X, theta_estimated = c(0, 0, 0),
    q_vec = q, T = 1L
  )
  expect_equal(fit[length(fit)], 2L)
})


test_that("public exact API returns NA for a singular full fit", {
  set.seed(92L)
  n <- 100L
  x <- rnorm(n)
  meta <- data.frame(x = x, duplicate_x = x)
  q <- rep(0.8, n)
  pic <- matrix(rbinom(n, size = 1L, prob = 0.4), nrow = 1L)
  rownames(pic) <- "peak1"

  expect_warning(
    fit <- pacs_test_cumu(
      covariate_meta.data = meta,
      formula_full = ~ x + duplicate_x,
      formula_null = ~ 1,
      pic_matrix = pic, max_T = 1L, cap_rates = q,
      n_cores = 1L, method = "exact"
    ),
    "non-converged"
  )
  expect_equal(unname(fit$pacs_converged[2L]), 2L)
  expect_true(is.na(fit$pacs_p_val))
})


test_that("exact optimizer distinguishes line-search and starting failures", {
  set.seed(1L)
  n <- 40L
  X <- matrix(rnorm(n), ncol = 1L)
  q <- runif(n, 0.3, 1)
  M <- as.integer(strsplit(
    "2002022200202021022010002000222000000012", "", fixed = TRUE
  )[[1L]])
  theta_start <- c(
    -3.12766210736320360,
     1.63829686962693089,
     0.77849518770127824
  )

  no_step <- PACS:::irls_iter_cumu(
    M_vec = M, X = X, theta_estimated = theta_start,
    q_vec = q, T = 2L, max_halving = 0L
  )
  expect_equal(no_step[length(no_step)], 4L)

  X_no_beta <- matrix(numeric(), nrow = 3L, ncol = 0L)
  invalid_start <- PACS:::irls_iter_cumu(
    M_vec = 0:2, X = X_no_beta,
    theta_estimated = c(40, -1000), q_vec = rep(1, 3L), T = 2L
  )
  expect_equal(invalid_start[length(invalid_start)], 5L)
})


test_that("converged exact fits also satisfy the score tolerance", {
  fx <- make_cumu_fixture(
    n = 250L, p = 2L, T = 2L,
    alpha = c(0.7, -0.3), beta = c(0.4, -0.2), seed = 811L
  )
  M <- PACS:::simulate_cumu_pacs(
    X = fx$X, alpha = fx$alpha, beta = fx$beta,
    q = fx$q, capture = "B", seed = 812L
  )
  theta_start <- PACS:::warm_start_theta(M, fx$q, T = 2L, p_beta = 2L)
  full <- PACS:::irls_iter_cumu(
    M, fx$X, theta_start, fx$q, T = 2L, score_tolerance = 1e-4
  )
  expect_equal(full[length(full)], 1L)
  expect_lte(max(abs(PACS:::loss_gradient_cumu(
    full[-length(full)], fx$X, M, fx$q, T = 2L
  ))), 1e-4)

  theta_null <- theta_start
  theta_null[4L] <- 0
  null <- PACS:::irls_iter_cumu_null(
    M, fx$X, theta_null, hold_zero = 4L,
    q_vec = fx$q, T = 2L, score_tolerance = 1e-4
  )
  expect_equal(null[length(null)], 1L)
  null_score <- PACS:::loss_gradient_cumu(
    null[-length(null)], fx$X, M, fx$q, T = 2L
  )
  expect_lte(max(abs(null_score[-4L])), 1e-4)
})


test_that("exact LRT withholds inference for failed fits", {
  n <- 20L
  X <- matrix(rep(c(-1, 1), n / 2L), ncol = 1L)
  M <- as.integer(X[, 1L] > 0)
  q <- rep(1, n)
  theta <- matrix(c(0, 0), ncol = 1L)

  expect_warning(
    p <- PACS:::compare_models_cumu(
      x_full = X, theta_estimated_full = theta,
      theta_estimated_null = theta, q_vec = q,
      c_by_r = matrix(M, ncol = 1L), T = 1L, df_test = 1L,
      conv_full = 2L, conv_null = 1L, mc.cores = 1L
    ),
    "non-converged"
  )
  expect_true(is.na(p))
})


test_that("materially negative exact LRT statistics produce NA", {
  n <- 20L
  X <- matrix(rep(c(-1, 1), n / 2L), ncol = 1L)
  M <- as.integer(X[, 1L] > 0)
  q <- rep(1, n)

  expect_warning(
    p <- PACS:::compare_models_cumu(
      x_full = X,
      theta_estimated_full = matrix(c(0, -5), ncol = 1L),
      theta_estimated_null = matrix(c(0, 0), ncol = 1L),
      q_vec = q, c_by_r = matrix(M, ncol = 1L),
      T = 1L, df_test = 1L, mc.cores = 1L
    ),
    "materially negative"
  )
  expect_true(is.na(p))
})


test_that("exact LRT checks both fits for order boundaries", {
  n <- 30L
  X <- matrix(rnorm(n), ncol = 1L)
  M <- rep(0:2, length.out = n)
  q <- rep(0.8, n)
  theta_full <- matrix(c(0, log(1), 0), ncol = 1L)
  theta_null <- matrix(c(0, log(1e-5), 0), ncol = 1L)

  expect_warning(
    p <- PACS:::compare_models_cumu(
      x_full = X, theta_estimated_full = theta_full,
      theta_estimated_null = theta_null, q_vec = q,
      c_by_r = matrix(M, ncol = 1L), T = 2L, df_test = 1L,
      mc.cores = 1L
    ),
    "order-constraint boundary"
  )
  expect_true(is.na(p))
})


test_that("exact LRT validates parameter dimensions", {
  X <- matrix(rnorm(20L), ncol = 1L)
  M <- matrix(rep(0:1, length.out = 20L), ncol = 1L)
  theta <- matrix(c(0, 0), ncol = 1L)

  expect_error(
    PACS:::compare_models_cumu(
      x_full = X, theta_estimated_full = theta[-1L, , drop = FALSE],
      theta_estimated_null = theta, q_vec = rep(0.8, 20L),
      c_by_r = M, T = 1L, df_test = 1L
    ),
    "must each have"
  )
})
