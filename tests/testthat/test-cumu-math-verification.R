## atilde_from_alpha_test() lives in helper-cumu.R.

finite_difference_score_test <- function(theta, X, M, q, T, h = 1e-6) {
  vapply(seq_along(theta), function(j) {
    theta_plus <- theta
    theta_minus <- theta
    theta_plus[j] <- theta_plus[j] + h
    theta_minus[j] <- theta_minus[j] - h
    (
      PACS:::loss_fun_cumu(theta_plus, X, M, q, T) -
        PACS:::loss_fun_cumu(theta_minus, X, M, q, T)
    ) / (2 * h)
  }, numeric(1L))
}


enumerate_information_test <- function(alpha, beta, X, q) {
  T <- length(alpha)
  n_parameters <- T + length(beta)
  information <- matrix(0, n_parameters, n_parameters)

  for (i in seq_len(nrow(X))) {
    X_i <- X[i, , drop = FALSE]
    link_i <- PACS:::cumu_link(alpha, beta, X_i)
    probabilities <- c(
      (1 - q[i]) + q[i] * link_i$one_minus_p[1L, 1L],
      q[i] * link_i$Delta[1L, ]
    )
    expect_equal(sum(probabilities), 1, tolerance = 1e-12)

    for (m in 0:T) {
      score_i <- PACS:::score_alpha_beta(
        alpha = alpha, beta = beta, X = X_i, M = m, q = q[i],
        link = link_i
      )
      information <- information +
        probabilities[m + 1L] * tcrossprod(score_i)
    }
  }
  information
}


test_that("exact score matches finite differences across numerical regimes", {
  set.seed(612L)
  cases <- list(
    interior = list(
      alpha = c(0.8, -0.1, -1.0), beta = c(0.4, -0.2),
      X = matrix(rnorm(36L), nrow = 18L, ncol = 2L),
      q = runif(18L, 0.35, 0.95), M = rep(0:3, length.out = 18L)
    ),
    sparse = list(
      alpha = c(-2.5, -4.0), beta = 0.3,
      X = matrix(rnorm(20L), ncol = 1L),
      q = runif(20L, 0.4, 0.9), M = c(rep(0L, 16L), 1L, 0L, 2L, 0L)
    ),
    no_capture_layer = list(
      alpha = c(0.5, -0.6), beta = -0.25,
      X = matrix(rnorm(15L), ncol = 1L),
      q = rep(1, 15L), M = rep(0:2, length.out = 15L)
    ),
    near_order_boundary = list(
      alpha = c(0.4, 0.399), beta = 0.1,
      X = matrix(rnorm(12L), ncol = 1L),
      q = runif(12L, 0.5, 1), M = rep(0:2, length.out = 12L)
    ),
    near_separation = list(
      alpha = c(40, 39.999), beta = numeric(),
      X = matrix(numeric(), nrow = 3L, ncol = 0L),
      q = rep(1, 3L), M = 0:2
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    T <- length(case$alpha)
    theta <- c(atilde_from_alpha_test(case$alpha), case$beta)
    analytic <- PACS:::loss_gradient_cumu(
      theta, case$X, case$M, case$q, T
    )
    numeric <- finite_difference_score_test(
      theta, case$X, case$M, case$q, T
    )
    relative_error <- max(
      abs(analytic - numeric) / pmax(1, abs(analytic), abs(numeric))
    )
    expect_lt(relative_error, 1e-6, label = case_name)
  }
})


test_that("expected information matches category enumeration", {
  set.seed(713L)
  cases <- list(
    interior = list(
      alpha = c(0.8, -0.1, -1.0), beta = c(0.4, -0.2),
      X = matrix(rnorm(24L), nrow = 12L, ncol = 2L),
      q = runif(12L, 0.35, 0.95)
    ),
    sparse = list(
      alpha = c(-2.5, -4.0), beta = 0.3,
      X = matrix(rnorm(12L), ncol = 1L),
      q = runif(12L, 0.4, 0.9)
    ),
    no_capture_layer = list(
      alpha = c(0.5, -0.6), beta = -0.25,
      X = matrix(rnorm(10L), ncol = 1L), q = rep(1, 10L)
    ),
    near_order_boundary = list(
      alpha = c(0.4, 0.399), beta = 0.1,
      X = matrix(rnorm(8L), ncol = 1L), q = runif(8L, 0.5, 1)
    ),
    near_separation = list(
      alpha = c(40, 39.999), beta = numeric(),
      X = matrix(numeric(), nrow = 3L, ncol = 0L), q = rep(1, 3L)
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    enumerated <- enumerate_information_test(
      case$alpha, case$beta, case$X, case$q
    )
    analytic <- PACS:::infor_mat_alpha_beta(
      case$alpha, case$beta, case$X, case$q
    )
    relative_error <- max(
      abs(analytic - enumerated) /
        pmax(1, abs(analytic), abs(enumerated))
    )
    expect_lt(relative_error, 1e-12, label = case_name)

    atilde <- atilde_from_alpha_test(case$alpha)
    theta <- c(atilde, case$beta)
    jacobian <- PACS:::full_jacobian(atilde, length(case$beta))
    transformed_enumeration <- crossprod(jacobian, enumerated) %*% jacobian
    transformed_analytic <- PACS:::infor_mat_cumu(
      theta, case$X, case$q, length(case$alpha)
    )
    expect_equal(
      transformed_analytic, transformed_enumeration,
      tolerance = 1e-11, info = case_name
    )
  }
})
