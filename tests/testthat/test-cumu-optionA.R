## Option A (fragment-level binomial thinning) for the exact cumulative-logit
## path. Math note notes/cumulative_logit_math.md §2, §3.2, §4.2, §5.2.
##
## The verification chain is deliberately layered so that no stage validates
## itself: an independently written brute-force mixture checks Pr(M = m);
## finite differences of that likelihood check the score; per-cell category
## enumeration of that score checks the expected information.


## --- independent reference implementations --------------------------------

## Pr(M_i = m) written straight from the note, in ordinary (non-log) space,
## sharing no code with the implementation.
reference_optA_prob <- function(alpha, beta, X, q) {
  T <- length(alpha)
  n <- nrow(X)
  xb <- if (ncol(X) > 0L) as.numeric(X %*% beta) else rep.int(0, n)
  cumulative <- matrix(0, nrow = n, ncol = T)
  for (t in seq_len(T)) {
    cumulative[, t] <- 1 / (1 + exp(-(alpha[t] + xb)))
  }
  category <- matrix(0, nrow = n, ncol = T + 1L)
  category[, 1L] <- 1 - cumulative[, 1L]
  if (T >= 2L) {
    for (k in seq_len(T - 1L)) {
      category[, k + 1L] <- cumulative[, k] - cumulative[, k + 1L]
    }
  }
  category[, T + 1L] <- cumulative[, T]

  out <- matrix(0, nrow = n, ncol = T + 1L)
  for (m in 0:T) {
    for (k in m:T) {
      out[, m + 1L] <- out[, m + 1L] +
        choose(k, m) * q^m * (1 - q)^(k - m) * category[, k + 1L]
    }
  }
  out
}


## Independently written Option A negative log-likelihood in the (atilde,
## beta) parameterisation, for use with a generic optimizer.
reference_optA_nll <- function(theta, X, M, q, T) {
  alpha <- theta[1L]
  if (T >= 2L) {
    alpha <- c(alpha, theta[1L] - cumsum(exp(theta[2:T])))
  }
  beta <- theta[(T + 1L):length(theta)]
  probabilities <- reference_optA_prob(alpha, beta, X, q)
  -sum(log(probabilities[cbind(seq_along(M), M + 1L)]))
}


optA_finite_difference_score <- function(theta, X, M, q, T, h = 1e-6) {
  vapply(seq_along(theta), function(j) {
    theta_plus <- theta
    theta_minus <- theta
    theta_plus[j] <- theta_plus[j] + h
    theta_minus[j] <- theta_minus[j] - h
    (
      PACS:::loss_fun_cumu(theta_plus, X, M, q, T, capture = "A") -
        PACS:::loss_fun_cumu(theta_minus, X, M, q, T, capture = "A")
    ) / (2 * h)
  }, numeric(1L))
}


## I(theta) = sum_i sum_m Pr(M_i = m) s_{i,m} s_{i,m}^T, one cell at a time.
optA_enumerate_information <- function(alpha, beta, X, q) {
  T <- length(alpha)
  n_parameters <- T + length(beta)
  information <- matrix(0, n_parameters, n_parameters)

  for (i in seq_len(nrow(X))) {
    X_i <- X[i, , drop = FALSE]
    link_i <- PACS:::cumu_link(alpha, beta, X_i)
    probabilities <- reference_optA_prob(alpha, beta, X_i, q[i])[1L, ]
    expect_equal(sum(probabilities), 1, tolerance = 1e-12)

    for (m in 0:T) {
      score_i <- PACS:::score_alpha_beta(
        alpha = alpha, beta = beta, X = X_i, M = m, q = q[i],
        link = link_i, capture = "A"
      )
      information <- information + probabilities[m + 1L] * tcrossprod(score_i)
    }
  }
  information
}


## Numerical regimes shared by the score and information tests. `sparse` and
## `near_order_boundary` mirror the Option B cases in
## test-cumu-math-verification.R; `unit_capture` pins the q = 1 edge where
## log(1 - q) is -Inf and the binomial weights must degenerate cleanly.
optA_cases <- function() {
  set.seed(4141L)
  list(
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
    unit_capture = list(
      alpha = c(0.5, -0.6), beta = -0.25,
      X = matrix(rnorm(15L), ncol = 1L),
      q = rep(1, 15L), M = rep(0:2, length.out = 15L)
    ),
    mixed_capture = list(
      alpha = c(0.5, -0.6), beta = -0.25,
      X = matrix(rnorm(15L), ncol = 1L),
      q = rep(c(1, 0.5, 0.05), length.out = 15L),
      M = rep(0:2, length.out = 15L)
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
    ),
    single_threshold = list(
      alpha = 0.4, beta = c(0.3, -0.2),
      X = matrix(rnorm(20L), nrow = 10L, ncol = 2L),
      q = runif(10L, 0.3, 1), M = rep(0:1, length.out = 10L)
    )
  )
}


## --- binomial thinning weights --------------------------------------------

test_that("thinning weights are the Binomial(k, q) pmf and handle q = 1", {
  q <- c(0.05, 0.4, 0.999, 1)
  T <- 3L
  log_w <- PACS:::thinning_log_weights(q, T)

  expect_equal(dim(log_w), c(length(q), T + 1L, T + 1L))
  expect_true(all(is.finite(log_w) | log_w == -Inf))
  expect_false(any(is.na(log_w)))

  for (i in seq_along(q)) {
    for (k in 0:T) {
      ## Column k of the array, over m = 0..k, must be dbinom(., k, q_i).
      weights <- exp(log_w[i, seq_len(k + 1L), k + 1L])
      expect_equal(weights, dbinom(0:k, size = k, prob = q[i]),
                   tolerance = 1e-12,
                   info = sprintf("q = %g, k = %d", q[i], k))
      expect_equal(sum(weights), 1, tolerance = 1e-12)
      ## Entries with m > k are impossible.
      if (k < T) {
        expect_true(all(log_w[i, (k + 2L):(T + 1L), k + 1L] == -Inf))
      }
    }
  }

  ## At q = 1 nothing is lost: the only surviving weight is m = k.
  expect_equal(exp(log_w[4L, , ]), diag(T + 1L), tolerance = 1e-12)
})


## --- likelihood ------------------------------------------------------------

test_that("Option A category probabilities match the brute-force mixture", {
  for (case_name in names(optA_cases())) {
    case <- optA_cases()[[case_name]]
    T <- length(case$alpha)
    link <- PACS:::cumu_link(case$alpha, case$beta, case$X)
    log_w <- PACS:::thinning_log_weights(case$q, T)
    analytic <- exp(PACS:::optA_log_prob(link, log_w, T))
    reference <- reference_optA_prob(case$alpha, case$beta, case$X, case$q)

    expect_equal(rowSums(analytic), rep.int(1, nrow(case$X)),
                 tolerance = 1e-12, info = case_name)
    expect_equal(analytic, reference, tolerance = 1e-10, info = case_name)
  }
})


test_that("Option A log-likelihood matches the brute-force mixture", {
  for (case_name in names(optA_cases())) {
    case <- optA_cases()[[case_name]]
    T <- length(case$alpha)
    theta <- c(atilde_from_alpha_test(case$alpha), case$beta)
    reference <- -reference_optA_nll(theta, case$X, case$M, case$q, T)
    analytic <- PACS:::loss_fun_cumu(theta, case$X, case$M, case$q, T,
                                     capture = "A")
    if (case_name == "near_separation") {
      ## The ordinary-space reference cannot be used here: at alpha_1 = 40,
      ## sigma(40) rounds to exactly 1, so its pi_{i0} underflows to 0 and the
      ## reference log-likelihood is -Inf. The implementation works in log
      ## space and must stay finite. Checked against a hand oracle below.
      expect_equal(reference, -Inf, info = case_name)
      expect_true(is.finite(analytic))
      next
    }
    expect_equal(analytic, reference, tolerance = 1e-10, info = case_name)
  }
})


test_that("Option A stays exact where the naive mixture underflows", {
  ## alpha = c(40, 39.999), q = 1, so Pr(M = m) = pi_m with no capture layer.
  ## Independent oracle built from log-scale logistic identities:
  ##   log pi_0 = log sigma(-alpha_1)
  ##   log pi_1 = log sigma(alpha_1) + log sigma(-alpha_2)
  ##              + log(1 - exp(alpha_2 - alpha_1))
  ##   log pi_2 = log sigma(alpha_2)
  alpha <- c(40, 39.999)
  X <- matrix(numeric(), nrow = 3L, ncol = 0L)
  q <- rep(1, 3L)
  oracle <- c(
    plogis(-alpha[1L], log.p = TRUE),
    plogis(alpha[1L], log.p = TRUE) + plogis(-alpha[2L], log.p = TRUE) +
      log1p(-exp(alpha[2L] - alpha[1L])),
    plogis(alpha[2L], log.p = TRUE)
  )
  expect_equal(sum(exp(oracle)), 1, tolerance = 1e-12)

  link <- PACS:::cumu_link(alpha, numeric(), X)
  log_prob <- PACS:::optA_log_prob(
    link, PACS:::thinning_log_weights(q, 2L), 2L
  )
  expect_equal(log_prob[1L, ], oracle, tolerance = 1e-12)

  ## Same quantity through the public likelihood entry point, one cell per
  ## observed category.
  theta <- c(atilde_from_alpha_test(alpha), numeric(0))
  for (m in 0:2) {
    expect_equal(
      PACS:::loss_fun_cumu(theta, X[1L, , drop = FALSE], m, 1, 2L,
                           capture = "A"),
      oracle[m + 1L], tolerance = 1e-12
    )
  }
})


## --- score -----------------------------------------------------------------

test_that("Option A score matches finite differences across regimes", {
  for (case_name in names(optA_cases())) {
    case <- optA_cases()[[case_name]]
    T <- length(case$alpha)
    theta <- c(atilde_from_alpha_test(case$alpha), case$beta)
    analytic <- PACS:::loss_gradient_cumu(
      theta, case$X, case$M, case$q, T, capture = "A"
    )
    numeric <- optA_finite_difference_score(
      theta, case$X, case$M, case$q, T
    )
    relative_error <- max(
      abs(analytic - numeric) / pmax(1, abs(analytic), abs(numeric))
    )
    expect_lt(relative_error, 1e-6, label = case_name)
  }
})


test_that("Option A score has mean zero under the fitted model", {
  for (case_name in names(optA_cases())) {
    case <- optA_cases()[[case_name]]
    T <- length(case$alpha)
    for (i in seq_len(nrow(case$X))) {
      X_i <- case$X[i, , drop = FALSE]
      link_i <- PACS:::cumu_link(case$alpha, case$beta, X_i)
      probabilities <- reference_optA_prob(
        case$alpha, case$beta, X_i, case$q[i]
      )[1L, ]
      expected_score <- rowSums(vapply(0:T, function(m) {
        probabilities[m + 1L] * PACS:::score_alpha_beta(
          alpha = case$alpha, beta = case$beta, X = X_i, M = m,
          q = case$q[i], link = link_i, capture = "A"
        )
      }, numeric(T + length(case$beta))))
      expect_lt(max(abs(expected_score)), 1e-10,
                label = sprintf("%s cell %d", case_name, i))
    }
  }
})


## --- expected information --------------------------------------------------

test_that("Option A expected information matches category enumeration", {
  for (case_name in names(optA_cases())) {
    case <- optA_cases()[[case_name]]
    enumerated <- optA_enumerate_information(
      case$alpha, case$beta, case$X, case$q
    )
    analytic <- PACS:::infor_mat_alpha_beta(
      case$alpha, case$beta, case$X, case$q, capture = "A"
    )
    relative_error <- max(
      abs(analytic - enumerated) / pmax(1, abs(analytic), abs(enumerated))
    )
    expect_lt(relative_error, 1e-12, label = case_name)
    expect_equal(analytic, t(analytic), tolerance = 1e-12, info = case_name)

    ## And the same after the order-map Jacobian transform.
    atilde <- atilde_from_alpha_test(case$alpha)
    theta <- c(atilde, case$beta)
    jacobian <- PACS:::full_jacobian(atilde, length(case$beta))
    expect_equal(
      PACS:::infor_mat_cumu(theta, case$X, case$q, length(case$alpha),
                            capture = "A"),
      crossprod(jacobian, enumerated) %*% jacobian,
      tolerance = 1e-11, info = case_name
    )
  }
})


test_that("Option A expected information equals minus the expected Hessian", {
  ## The information equality I(theta) = -E[d^2 ell / d theta d theta'] holds
  ## only if the score and the category probabilities are mutually consistent,
  ## so this catches errors that the score-squared enumeration cannot: the
  ## enumeration would reproduce a wrong score's outer product happily.
  set.seed(6201L)
  alpha <- c(0.9, -0.2, -1.1)
  beta <- c(0.45, -0.3)
  X <- matrix(rnorm(8L), nrow = 4L, ncol = 2L)
  q <- c(0.4, 0.65, 0.9, 1)
  T <- length(alpha)
  h <- 1e-5

  negative_expected_hessian <- matrix(0, T + length(beta), T + length(beta))
  for (i in seq_len(nrow(X))) {
    X_i <- X[i, , drop = FALSE]
    probabilities <- reference_optA_prob(alpha, beta, X_i, q[i])[1L, ]
    for (m in 0:T) {
      ## Central-difference Jacobian of the score in (alpha, beta).
      hessian <- vapply(seq_len(T + length(beta)), function(j) {
        bump <- function(sign) {
          parameters <- c(alpha, beta)
          parameters[j] <- parameters[j] + sign * h
          PACS:::score_alpha_beta(
            alpha = parameters[seq_len(T)],
            beta = parameters[(T + 1L):length(parameters)],
            X = X_i, M = m, q = q[i], capture = "A"
          )
        }
        (bump(1) - bump(-1)) / (2 * h)
      }, numeric(T + length(beta)))
      negative_expected_hessian <- negative_expected_hessian -
        probabilities[m + 1L] * hessian
    }
  }

  analytic <- PACS:::infor_mat_alpha_beta(alpha, beta, X, q, capture = "A")
  relative_error <- max(
    abs(analytic - negative_expected_hessian) /
      pmax(1e-3, abs(analytic), abs(negative_expected_hessian))
  )
  expect_lt(relative_error, 1e-5)
})


test_that("Option A expected information is positive semi-definite", {
  for (case_name in names(optA_cases())) {
    case <- optA_cases()[[case_name]]
    information <- PACS:::infor_mat_alpha_beta(
      case$alpha, case$beta, case$X, case$q, capture = "A"
    )
    eigenvalues <- eigen(
      (information + t(information)) / 2, symmetric = TRUE, only.values = TRUE
    )$values
    expect_gt(min(eigenvalues), -1e-10 * max(1, max(eigenvalues)))
  }
})


test_that("Option A information is not tridiagonal in alpha when T >= 3", {
  ## Under Option B a cell contributes to alpha_t only for m in {t-1, t}, so
  ## the alpha block is tridiagonal. Option A spreads posterior mass over all
  ## latent k >= m, which couples non-adjacent thresholds. This guards the
  ## note's §5.2 claim that the closed-form Option B blocks cannot be reused.
  case <- optA_cases()$interior
  information <- PACS:::infor_mat_alpha_beta(
    case$alpha, case$beta, case$X, case$q, capture = "A"
  )
  expect_gt(abs(information[1L, 3L]), 1e-6)

  information_B <- PACS:::infor_mat_alpha_beta(
    case$alpha, case$beta, case$X, case$q, capture = "B"
  )
  expect_equal(information_B[1L, 3L], 0, tolerance = 1e-14)
})


## --- reduction to Option B and to the binary model at T = 1 ---------------

test_that("Option A and Option B coincide at T = 1", {
  fx <- make_cumu_fixture(n = 120L, p = 2L, T = 1L,
                          alpha = 0.45, beta = c(0.5, -0.3), seed = 8101L)
  ## Include an exactly-complete-capture cell.
  fx$q[1L] <- 1
  M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                 q = fx$q, capture = "A", seed = 8102L)
  expect_true(all(M %in% 0:1))
  theta <- c(fx$alpha, fx$beta)

  ll_A <- PACS:::loss_fun_cumu(theta, fx$X, M, fx$q, 1L, capture = "A")
  ll_B <- PACS:::loss_fun_cumu(theta, fx$X, M, fx$q, 1L, capture = "B")
  ## Option B drops the constant log q_i for m_i >= 1; Option A keeps it.
  expect_equal(ll_A, ll_B + sum(log(fx$q[M == 1L])), tolerance = 1e-10)

  expect_equal(
    PACS:::loss_gradient_cumu(theta, fx$X, M, fx$q, 1L, capture = "A"),
    PACS:::loss_gradient_cumu(theta, fx$X, M, fx$q, 1L, capture = "B"),
    tolerance = 1e-10
  )
  expect_equal(
    PACS:::infor_mat_cumu(theta, fx$X, fx$q, 1L, capture = "A"),
    PACS:::infor_mat_cumu(theta, fx$X, fx$q, 1L, capture = "B"),
    tolerance = 1e-10
  )
})


test_that("Option A reduces to the binary PACS model at T = 1", {
  fx <- make_cumu_fixture(n = 120L, p = 2L, T = 1L,
                          alpha = 0.45, beta = c(0.5, -0.3), seed = 8201L)
  fx$q[1L] <- 1
  M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                 q = fx$q, capture = "A", seed = 8202L)
  theta <- c(fx$alpha, fx$beta)
  X_binary <- cbind(intercept = 1, fx$X)
  p_bg <- plogis(as.numeric(X_binary %*% theta))

  ## Option A keeps every capture term, so it equals the binary likelihood
  ## exactly rather than up to a constant.
  expect_equal(
    PACS:::loss_fun_cumu(theta, fx$X, M, fx$q, 1L, capture = "A"),
    loss_fun(p_bg = p_bg, q_vec = fx$q, y_vec = M),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(PACS:::loss_gradient_cumu(theta, fx$X, M, fx$q, 1L,
                                         capture = "A")),
    as.numeric(loss_gradient(xdumm = X_binary, p_bg = p_bg, q_vec = fx$q,
                             y_vec = M)),
    tolerance = 1e-10
  )
  expect_equal(
    as.numeric(PACS:::infor_mat_cumu(theta, fx$X, fx$q, 1L, capture = "A")),
    as.numeric(infor_mat(xdumm = X_binary, p_bg = p_bg, q_vec = fx$q)),
    tolerance = 1e-10
  )
})


test_that("Option A and Option B coincide at q = 1 for any T", {
  ## With no capture layer both models reduce to the complete-data cumulative
  ## logit, so they must agree exactly at every T -- not just at T = 1. This
  ## also restores tridiagonality of the alpha block, since the thinning
  ## weights collapse to 1[m = k].
  set.seed(8001L)
  for (T in 2:4) {
    n <- 12L
    X <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)
    q <- rep(1, n)
    alpha <- sort(rnorm(T), decreasing = TRUE)
    beta <- rnorm(2L)
    theta <- c(atilde_from_alpha_test(alpha), beta)
    M <- rep(0:T, length.out = n)

    expect_equal(
      PACS:::loss_fun_cumu(theta, X, M, q, T, capture = "A"),
      PACS:::loss_fun_cumu(theta, X, M, q, T, capture = "B"),
      tolerance = 1e-12, info = sprintf("T = %d", T)
    )
    expect_equal(
      PACS:::loss_gradient_cumu(theta, X, M, q, T, capture = "A"),
      PACS:::loss_gradient_cumu(theta, X, M, q, T, capture = "B"),
      tolerance = 1e-12, info = sprintf("T = %d", T)
    )
    information <- PACS:::infor_mat_alpha_beta(alpha, beta, X, q,
                                               capture = "A")
    expect_equal(
      PACS:::infor_mat_cumu(theta, X, q, T, capture = "A"),
      PACS:::infor_mat_cumu(theta, X, q, T, capture = "B"),
      tolerance = 1e-12, info = sprintf("T = %d", T)
    )
    if (T >= 3L) {
      ## Tridiagonal: thresholds more than one apart do not couple.
      expect_equal(max(abs(information[1L, 3:T])), 0,
                   tolerance = 1e-14, info = sprintf("T = %d", T))
    }
  }
})


## --- simulator -------------------------------------------------------------

test_that("Option A simulator matches the analytic observed distribution", {
  set.seed(8301L)
  n <- 40000L
  alpha <- c(0.7, -0.5)
  beta <- 0.6
  X <- matrix(rep(c(-1, 1), length.out = n), ncol = 1L)
  q <- rep(c(0.35, 0.85), length.out = n)

  M <- PACS:::simulate_cumu_pacs(X = X, alpha = alpha, beta = beta, q = q,
                                 capture = "A", seed = 8302L)
  expect_true(all(M %in% 0:2))

  analytic <- reference_optA_prob(alpha, beta, X, q)
  ## Two distinct (x, q) cells. Within each stratum the counts are multinomial
  ## with the analytic probabilities, so a goodness-of-fit test is the honest
  ## comparison — an absolute tolerance would either be flaky on the rare
  ## M = 2 category or too loose to detect a wrong simulator.
  for (stratum in 1:2) {
    rows <- seq(stratum, n, by = 2L)
    observed <- as.numeric(table(factor(M[rows], levels = 0:2)))
    expected <- analytic[stratum, ] * length(rows)
    statistic <- sum((observed - expected)^2 / expected)
    expect_gt(pchisq(statistic, df = 2L, lower.tail = FALSE), 1e-3)
  }

  ## Under Option A a partially captured cell can still report a nonzero
  ## count; under Option B its count collapses to zero. So Option A puts less
  ## mass on M = 0 for the same (alpha, beta, q).
  M_B <- PACS:::simulate_cumu_pacs(X = X, alpha = alpha, beta = beta, q = q,
                                   capture = "B", seed = 8302L)
  expect_lt(mean(M == 0L), mean(M_B == 0L))
})


## --- estimation ------------------------------------------------------------

test_that("Option A Fisher scoring reaches the same optimum as a generic optimizer", {
  skip_on_cran()
  fx <- make_cumu_fixture(n = 600L, p = 2L, T = 2L,
                          alpha = c(0.8, -0.6), beta = c(0.6, -0.35),
                          seed = 8401L)
  M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                 q = fx$q, capture = "A", seed = 8402L)
  theta_init <- PACS:::warm_start_theta(M = M, q = fx$q, T = 2L, p_beta = 2L)

  fit <- PACS:::irls_iter_cumu(M_vec = M, X = fx$X,
                               theta_estimated = theta_init,
                               q_vec = fx$q, T = 2L, capture = "A")
  expect_equal(fit[length(fit)], 1L)
  theta_hat <- fit[seq_len(length(fit) - 1L)]

  reference_fit <- optim(
    par = theta_init, fn = reference_optA_nll,
    X = fx$X, M = M, q = fx$q, T = 2L,
    method = "BFGS", control = list(reltol = 1e-14, maxit = 2000L)
  )
  expect_equal(reference_fit$convergence, 0L)
  expect_equal(theta_hat, reference_fit$par, tolerance = 1e-4)
  expect_lte(
    max(abs(PACS:::loss_gradient_cumu(theta_hat, fx$X, M, fx$q, 2L,
                                      capture = "A"))),
    1e-4
  )
})


test_that("Option A recovers alpha and beta on Option A data", {
  skip_on_cran()
  n_rep <- 30L
  alpha <- c(0.8, -0.6)
  beta <- c(0.6, -0.35)

  alpha_hat <- matrix(NA_real_, nrow = n_rep, ncol = 2L)
  beta_hat <- matrix(NA_real_, nrow = n_rep, ncol = 2L)
  conv <- integer(n_rep)

  for (r in seq_len(n_rep)) {
    fx <- make_cumu_fixture(n = 400L, p = 2L, T = 2L, alpha = alpha,
                            beta = beta, seed = 8500L + r)
    M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                   q = fx$q, capture = "A", seed = 8600L + r)
    theta_init <- PACS:::warm_start_theta(M = M, q = fx$q, T = 2L,
                                          p_beta = 2L)
    res <- PACS:::irls_iter_cumu(M_vec = M, X = fx$X,
                                 theta_estimated = theta_init,
                                 q_vec = fx$q, T = 2L, capture = "A")
    conv[r] <- res[length(res)]
    if (conv[r] == 1L) {
      th <- res[seq_len(length(res) - 1L)]
      alpha_hat[r, ] <- PACS:::alpha_from_atilde(th[1:2])
      beta_hat[r, ] <- th[3:4]
    }
  }

  ok <- conv == 1L
  expect_gt(mean(ok), 0.85)
  expect_lt(max(abs(colMeans(alpha_hat[ok, , drop = FALSE]) - alpha)), 0.15)
  expect_lt(max(abs(colMeans(beta_hat[ok, , drop = FALSE]) - beta)), 0.15)
})


test_that("the capture model is not interchangeable at T = 2", {
  skip_on_cran()
  ## Fitting Option B to fragment-thinned data (and vice versa) distorts the
  ## thresholds systematically, which is the whole reason Option A exists.
  ## The note (§2) predicts the two agree only when p_{i2} is negligible.
  n_rep <- 25L
  alpha <- c(0.8, -0.6)
  beta <- c(0.6, -0.35)

  alpha_matched <- matrix(NA_real_, nrow = n_rep, ncol = 2L)
  alpha_mismatched <- matrix(NA_real_, nrow = n_rep, ncol = 2L)

  for (r in seq_len(n_rep)) {
    fx <- make_cumu_fixture(n = 500L, p = 2L, T = 2L, alpha = alpha,
                            beta = beta, seed = 8700L + r)
    M <- PACS:::simulate_cumu_pacs(X = fx$X, alpha = fx$alpha, beta = fx$beta,
                                   q = fx$q, capture = "A", seed = 8800L + r)
    theta_init <- PACS:::warm_start_theta(M = M, q = fx$q, T = 2L,
                                          p_beta = 2L)
    matched <- PACS:::irls_iter_cumu(M_vec = M, X = fx$X,
                                     theta_estimated = theta_init,
                                     q_vec = fx$q, T = 2L, capture = "A")
    mismatched <- PACS:::irls_iter_cumu(M_vec = M, X = fx$X,
                                        theta_estimated = theta_init,
                                        q_vec = fx$q, T = 2L, capture = "B")
    if (matched[length(matched)] == 1L) {
      alpha_matched[r, ] <- PACS:::alpha_from_atilde(matched[1:2])
    }
    if (mismatched[length(mismatched)] == 1L) {
      alpha_mismatched[r, ] <- PACS:::alpha_from_atilde(mismatched[1:2])
    }
  }

  ok <- stats::complete.cases(alpha_matched) &
    stats::complete.cases(alpha_mismatched)
  expect_gt(sum(ok), 15L)

  bias_matched <- colMeans(alpha_matched[ok, , drop = FALSE]) - alpha
  bias_mismatched <- colMeans(alpha_mismatched[ok, , drop = FALSE]) - alpha
  expect_lt(max(abs(bias_matched)), 0.15)
  expect_gt(max(abs(bias_mismatched)), 0.30)
})


## --- public API ------------------------------------------------------------

test_that("pacs_test_cumu exposes Option A and defaults to Option B", {
  set.seed(8901L)
  n <- 300L
  n_peaks <- 5L
  group <- factor(rep(c("ctrl", "treat"), each = n / 2L))
  meta <- data.frame(group = group)
  X_raw <- model.matrix(~group, data = meta)[, -1, drop = FALSE]
  q <- runif(n, 0.35, 0.9)

  pic <- matrix(0L, nrow = n_peaks, ncol = n)
  rownames(pic) <- paste0("peak", seq_len(n_peaks))
  for (j in seq_len(n_peaks)) {
    pic[j, ] <- PACS:::simulate_cumu_pacs(
      X = X_raw, alpha = c(0.7, -0.5), beta = 0.8, q = q,
      capture = "A", seed = 8900L + j
    )
  }

  result_A <- pacs_test_cumu(
    covariate_meta.data = meta, formula_full = ~group, formula_null = ~1,
    pic_matrix = pic, max_T = 2, cap_rates = q, n_cores = 1L,
    method = "exact", capture = "A"
  )
  expect_named(result_A, c("pacs_converged", "pacs_p_val"))
  expect_length(result_A$pacs_p_val, n_peaks)
  expect_equal(names(result_A$pacs_p_val), rownames(pic))
  expect_true(all(result_A$pacs_p_val >= 0 & result_A$pacs_p_val <= 1,
                  na.rm = TRUE))
  expect_true(all(result_A$pacs_converged == 1L))

  ## The default must remain Option B, unchanged by this feature.
  result_default <- pacs_test_cumu(
    covariate_meta.data = meta, formula_full = ~group, formula_null = ~1,
    pic_matrix = pic, max_T = 2, cap_rates = q, n_cores = 1L,
    method = "exact"
  )
  result_B <- pacs_test_cumu(
    covariate_meta.data = meta, formula_full = ~group, formula_null = ~1,
    pic_matrix = pic, max_T = 2, cap_rates = q, n_cores = 1L,
    method = "exact", capture = "B"
  )
  expect_equal(result_default, result_B)
  ## And Option A must actually change the answer on thinned data.
  expect_false(isTRUE(all.equal(result_A$pacs_p_val, result_B$pacs_p_val)))

  expect_error(
    pacs_test_cumu(
      covariate_meta.data = meta, formula_full = ~group, formula_null = ~1,
      pic_matrix = pic, max_T = 2, cap_rates = q, n_cores = 1L,
      method = "exact", capture = "fragment"
    ),
    "should be one of"
  )

  ## The stacked path has no capture model; asking for one must not pass
  ## unnoticed.
  expect_warning(
    pacs_test_cumu(
      covariate_meta.data = meta, formula_full = ~group, formula_null = ~1,
      pic_matrix = pic, max_T = 2, cap_rates = q, n_cores = 1L,
      method = "stacked", capture = "A"
    ),
    "only applies to method"
  )
})


test_that("Option A warns when it top-codes counts, Option B does not", {
  ## Capping the observed count commutes with cell-level dropout but not with
  ## fragment thinning, so the same top-code is exact under B and an
  ## approximation under A. Other warnings (convergence, boundary) are
  ## unrelated, so match on the top-coding message specifically.
  top_code_warnings <- function(...) {
    messages <- character(0)
    withCallingHandlers(
      pacs_test_cumu(...),
      warning = function(w) {
        messages <<- c(messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    grep("exceed max_T", messages, value = TRUE)
  }

  set.seed(9301L)
  n <- 300L
  n_peaks <- 3L
  group <- factor(rep(c("ctrl", "treat"), each = n / 2L))
  meta <- data.frame(group = group)
  X_raw <- model.matrix(~group, data = meta)[, -1, drop = FALSE]
  q <- runif(n, 0.5, 0.95)

  ## Simulate with three latent thresholds, then fit with max_T = 2 so that
  ## genuine counts of 3 must be top-coded.
  pic <- matrix(0L, nrow = n_peaks, ncol = n)
  rownames(pic) <- paste0("peak", seq_len(n_peaks))
  for (j in seq_len(n_peaks)) {
    pic[j, ] <- PACS:::simulate_cumu_pacs(
      X = X_raw, alpha = c(1.4, 0.5, -0.4), beta = 0.6, q = q,
      capture = "A", seed = 9310L + j
    )
  }
  expect_gt(sum(pic > 2L), 0L)

  expect_match(
    top_code_warnings(
      covariate_meta.data = meta, formula_full = ~group, formula_null = ~1,
      pic_matrix = pic, max_T = 2, cap_rates = q, n_cores = 1L,
      method = "exact", capture = "A"
    ),
    sprintf("%d observed count", sum(pic > 2L))
  )
  ## Under Option B the top-code is exact, so it is not flagged.
  expect_length(
    top_code_warnings(
      covariate_meta.data = meta, formula_full = ~group, formula_null = ~1,
      pic_matrix = pic, max_T = 2, cap_rates = q, n_cores = 1L,
      method = "exact", capture = "B"
    ),
    0L
  )
  ## And nothing is flagged under A when max_T already covers the data.
  expect_length(
    top_code_warnings(
      covariate_meta.data = meta, formula_full = ~group, formula_null = ~1,
      pic_matrix = pic, max_T = 3, cap_rates = q, n_cores = 1L,
      method = "exact", capture = "A"
    ),
    0L
  )
})


test_that("Option A LRT is calibrated under the null", {
  skip_on_cran()
  n_rep <- 200L
  n <- 400L
  T <- 2L
  alpha <- c(0.7, -0.4)

  p_values <- rep(NA_real_, n_rep)
  for (r in seq_len(n_rep)) {
    set.seed(9100L + r)
    X <- matrix(rnorm(n * 2L), nrow = n, ncol = 2L)
    colnames(X) <- c("x1", "x2")
    q <- runif(n, 0.4, 0.9)
    ## Second covariate has no effect; it is the one tested.
    M <- PACS:::simulate_cumu_pacs(X = X, alpha = alpha, beta = c(0.4, 0),
                                   q = q, capture = "A")

    theta_init <- PACS:::warm_start_theta(M = M, q = q, T = T, p_beta = 2L)
    theta_init_null <- theta_init
    theta_init_null[T + 2L] <- 0
    full_fit <- PACS:::irls_iter_cumu(
      M_vec = M, X = X, theta_estimated = theta_init, q_vec = q, T = T,
      capture = "A"
    )
    null_fit <- PACS:::irls_iter_cumu_null(
      M_vec = M, X = X, theta_estimated = theta_init_null,
      hold_zero = T + 2L, q_vec = q, T = T, capture = "A"
    )
    if (full_fit[length(full_fit)] != 1L ||
        null_fit[length(null_fit)] != 1L) {
      next
    }
    p_values[r] <- PACS:::compare_models_cumu(
      x_full = X,
      theta_estimated_full = matrix(full_fit[1:4], ncol = 1L),
      theta_estimated_null = matrix(null_fit[1:4], ncol = 1L),
      q_vec = q, c_by_r = matrix(M, ncol = 1L), T = T, df_test = 1L,
      capture = "A", mc.cores = 1L
    )
  }

  usable <- !is.na(p_values)
  expect_gt(mean(usable), 0.9)
  rejection_rate <- mean(p_values[usable] < 0.05)
  ## Binomial SE at 200 replicates is about 0.015; allow ~3 SE.
  expect_lt(rejection_rate, 0.10)
  expect_gt(mean(p_values[usable] < 0.5), 0.35)
})
