## Cumulative-logit PACS, exact (Option B)
##
## Implements the proper proportional-odds likelihood with cell-level
## all-or-nothing capture (notes/cumulative_logit_math.md, §3.1, §4.1, §5.1,
## §6, §7). Internal helpers; the public surface is pacs_test_cumu.
##
## At T = 1 the likelihood, score, and expected information reduce
## term-for-term to loss_fun, loss_gradient, and infor_mat in
## R/param-estimate_logit_get_p_by_t_June.R; tests in
## tests/testthat/test-cumu-T1-parity.R verify this.


is_error_cumu <- function(x) inherits(x, "try-error")


validate_cumu_response <- function(M, T) {
  if (length(T) != 1L || !is.finite(T) || T < 1L || T != as.integer(T)) {
    stop("T must be a positive integer.", call. = FALSE)
  }
  if (any(!is.finite(M)) || any(M < 0) || any(M != floor(M))) {
    stop("M must contain finite, non-negative integer counts.", call. = FALSE)
  }
  if (any(M > T)) {
    stop(sprintf(
      "M must be in 0:%d. Top-code counts above T before calling this internal function.",
      as.integer(T)
    ), call. = FALSE)
  }
  invisible(TRUE)
}


## log(exp(a) + exp(b)), including the case where either input is -Inf.
logspace_add_cumu <- function(a, b) {
  hi <- pmax(a, b)
  lo <- pmin(a, b)
  out <- hi + log1p(exp(lo - hi))
  both_neg_inf <- is.infinite(hi) & hi < 0
  out[both_neg_inf] <- -Inf
  out
}


## log(1 - exp(x)) for x <= 0, stable when x is close to zero.
log1mexp_cumu <- function(x) {
  out <- numeric(length(x))
  near_zero <- x > -log(2)
  out[near_zero] <- log(-expm1(x[near_zero]))
  out[!near_zero] <- log1p(-exp(x[!near_zero]))
  out
}


## --- reparameterisation ---------------------------------------------------

## ã -> α with α_1 = ã_1, α_t = ã_1 - sum_{s=2..t} exp(ã_s).
## Strictly enforces α_1 >= α_2 >= ... >= α_T.
alpha_from_atilde <- function(atilde) {
  T <- length(atilde)
  if (T == 1L) return(atilde)
  alpha <- numeric(T)
  alpha[1] <- atilde[1]
  cs <- cumsum(exp(atilde[2:T]))
  alpha[2:T] <- atilde[1] - cs
  alpha
}

## Jacobian J = ∂α/∂ã. Lower-triangular: J[,1] = 1, J[t,s] = -exp(ã_s) for
## 2 <= s <= t, zero otherwise.
atilde_jacobian <- function(atilde) {
  T <- length(atilde)
  J <- matrix(0, nrow = T, ncol = T)
  J[, 1] <- 1
  if (T >= 2L) {
    e <- exp(atilde[2:T])
    for (s in 2:T) {
      J[s:T, s] <- -e[s - 1]
    }
  }
  J
}

## Block Jacobian for θ = (ã, β): blkdiag(J, I_p).
full_jacobian <- function(atilde, p_beta) {
  T <- length(atilde)
  Jt <- atilde_jacobian(atilde)
  m <- T + p_beta
  J <- matrix(0, nrow = m, ncol = m)
  J[1:T, 1:T] <- Jt
  if (p_beta > 0L) {
    diag(J)[(T + 1):m] <- 1
  }
  J
}


## --- link-function quantities --------------------------------------------

## Compute p_{it}, u_{it}, Delta_{it}, D_{it} given current (α, β).
## Conventions: p_{i,0} = 1, p_{i,T+1} = 0, u_{i,0} = u_{i,T+1} = 0.
##
## Returns a list with n×T matrices p, u, Delta, D and their stable log/ratio
## counterparts. Delta is evaluated without subtracting two rounded logistic
## probabilities.
cumu_link <- function(alpha, beta, X) {
  n <- nrow(X)
  T <- length(alpha)
  if (length(beta) > 0L) {
    xb <- as.numeric(X %*% beta)
  } else {
    xb <- rep.int(0, n)
  }
  ## eta[i, t] = α_t + x_i^T β, so p[i, t] = σ(eta).
  eta <- outer(xb, alpha, FUN = "+")
  p <- plogis(eta)
  one_minus_p <- plogis(eta, lower.tail = FALSE)
  log_p <- plogis(eta, log.p = TRUE)
  log_one_minus_p <- plogis(eta, lower.tail = FALSE, log.p = TRUE)
  log_u <- log_p + log_one_minus_p
  u <- exp(log_u)

  ## For a >= b,
  ##   log{sigma(a) - sigma(b)}
  ##     = log sigma(a) + log sigma(-b) + log{1 - exp(b-a)}.
  ## Here a-b is the threshold gap and does not depend on x_i^T beta.
  log_Delta <- matrix(0, nrow = n, ncol = T)
  if (T >= 2L) {
    for (t in seq_len(T - 1L)) {
      gap <- alpha[t] - alpha[t + 1L]
      log_Delta[, t] <- log_p[, t] + log_one_minus_p[, t + 1L] +
        log1mexp_cumu(-gap)
    }
    log_Delta[, T] <- log_p[, T]
  } else {
    log_Delta[, 1L] <- log_p[, 1L]
  }
  Delta <- exp(log_Delta)

  ## D/Delta = 1 - p_t - p_{t+1}; evaluate this ratio directly so the score
  ## remains finite even when both logistic probabilities round to one.
  D_over_Delta <- matrix(0, nrow = n, ncol = T)
  if (T >= 2L) {
    D_over_Delta[, seq_len(T - 1L)] <-
      one_minus_p[, seq_len(T - 1L), drop = FALSE] -
      p[, 2:T, drop = FALSE]
  }
  D_over_Delta[, T] <- one_minus_p[, T]
  D <- Delta * D_over_Delta

  list(
    p = p, one_minus_p = one_minus_p,
    log_p = log_p, log_one_minus_p = log_one_minus_p,
    u = u, log_u = log_u,
    Delta = Delta, log_Delta = log_Delta,
    D = D, D_over_Delta = D_over_Delta
  )
}


## --- log-likelihood (Option B) -------------------------------------------

## Per-cell log Pr(M_i | x_i, q_i) under Option B (note §3.1). Drops the
## constant log q_i term for m_i >= 1 (irrelevant to optimisation).
loss_fun_cumu <- function(theta, X, M, q, T) {
  validate_cumu_response(M, T)
  alpha <- alpha_from_atilde(theta[1:T])
  beta <- if (length(theta) > T) theta[(T + 1):length(theta)] else numeric(0)
  L <- cumu_link(alpha, beta, X)
  ll <- numeric(length(M))
  zero <- (M == 0L)
  if (any(zero)) {
    log_prob_zero <- logspace_add_cumu(
      log1p(-q[zero]),
      log(q[zero]) + L$log_one_minus_p[zero, 1L]
    )
    ll[zero] <- log_prob_zero
  }
  if (any(!zero)) {
    idx <- which(!zero)
    ll[idx] <- L$log_Delta[cbind(idx, M[idx])]
  }
  sum(ll)
}


## --- score in (α, β) -----------------------------------------------------

## Returns a length-(T + p_beta) vector: ∂ℓ / ∂(α, β) evaluated at
## (alpha, beta). Used internally; loss_gradient_cumu wraps with the
## reparameterisation.
score_alpha_beta <- function(alpha, beta, X, M, q, link = NULL) {
  T <- length(alpha)
  p_beta <- length(beta)
  if (is.null(link)) link <- cumu_link(alpha, beta, X)
  validate_cumu_response(M, T)
  u <- link$u

  s_alpha <- numeric(T)
  s_beta <- numeric(p_beta)

  zero <- (M == 0L)
  ## ----- m_i = 0 -----
  if (any(zero)) {
    iz <- which(zero)
    log_prob_zero <- logspace_add_cumu(
      log1p(-q[iz]),
      log(q[iz]) + link$log_one_minus_p[iz, 1L]
    )
    coef_zero <- -exp(log(q[iz]) + link$log_u[iz, 1L] - log_prob_zero)
    s_alpha[1] <- s_alpha[1] + sum(coef_zero)
    if (p_beta > 0L) {
      s_beta <- s_beta + as.numeric(crossprod(X[iz, , drop = FALSE], coef_zero))
    }
  }
  ## ----- m_i >= 1 -----
  if (any(!zero)) {
    inz <- which(!zero)
    Mi <- M[inz]
    ## ∂ℓ/∂β = (D_{i,M_i} / Delta_{i,M_i}) x_i
    coef_beta <- link$D_over_Delta[cbind(inz, Mi)]
    if (p_beta > 0L) {
      s_beta <- s_beta + as.numeric(crossprod(X[inz, , drop = FALSE], coef_beta))
    }
    ## ∂ℓ/∂α_t = u_{it}/Delta_{i,M_i} · 1[t=M_i] - u_{it}/Delta_{i,M_i} · 1[t=M_i+1]
    ## Group cells by Mi for vectorised accumulation.
    for (m in unique(Mi)) {
      sub <- which(Mi == m)
      cells <- inz[sub]
      ## t = m contribution (positive)
      ratio_left <- exp(
        link$log_u[cells, m] - link$log_Delta[cells, m]
      )
      s_alpha[m] <- s_alpha[m] + sum(ratio_left)
      ## t = m + 1 contribution (negative), only if m + 1 <= T
      if ((m + 1L) <= T) {
        ratio_right <- exp(
          link$log_u[cells, m + 1L] - link$log_Delta[cells, m]
        )
        s_alpha[m + 1L] <- s_alpha[m + 1L] - sum(ratio_right)
      }
    }
  }

  c(s_alpha, s_beta)
}


## Score in θ = (ã, β): apply the Jacobian transform score_θ = J̃^T · score_(α,β).
loss_gradient_cumu <- function(theta, X, M, q, T) {
  atilde <- theta[1:T]
  alpha <- alpha_from_atilde(atilde)
  beta <- if (length(theta) > T) theta[(T + 1):length(theta)] else numeric(0)
  s_ab <- score_alpha_beta(alpha, beta, X, M, q)
  J <- full_jacobian(atilde, length(beta))
  as.numeric(crossprod(J, s_ab))
}


## --- Fisher information in (α, β) (Option B, §5.1) -----------------------

infor_mat_alpha_beta <- function(alpha, beta, X, q, link = NULL) {
  T <- length(alpha)
  p_beta <- length(beta)
  if (is.null(link)) link <- cumu_link(alpha, beta, X)
  u <- link$u

  ## ν_{i,0} = q^2 u_{i,1}^2 / (1 - q p_{i,1}); contribution from m=0.
  log_prob_zero <- logspace_add_cumu(
    log1p(-q), log(q) + link$log_one_minus_p[, 1L]
  )
  nu0 <- exp(2 * log(q) + 2 * link$log_u[, 1L] - log_prob_zero)

  ## Stable forms of u_t^2 / Delta_t and u_t u_{t+1} / Delta_t.
  u2_over_delta <- exp(2 * link$log_u - link$log_Delta)

  ## ----- α-block (tridiagonal) -----
  I_aa <- matrix(0, nrow = T, ncol = T)
  ## diagonal
  I_aa[1, 1] <- sum(q * u2_over_delta[, 1]) + sum(nu0)
  if (T >= 2L) {
    for (t in 2:T) {
      right_term <- exp(2 * link$log_u[, t] - link$log_Delta[, t - 1L])
      I_aa[t, t] <- sum(q * u2_over_delta[, t]) + sum(q * right_term)
    }
    ## off-diagonal (super-diagonal; symmetric copy below)
    for (t in 1:(T - 1L)) {
      cross_term <- exp(
        link$log_u[, t] + link$log_u[, t + 1L] - link$log_Delta[, t]
      )
      val <- -sum(q * cross_term)
      I_aa[t, t + 1L] <- val
      I_aa[t + 1L, t] <- val
    }
  }

  ## ----- β-block -----
  if (p_beta > 0L) {
    v <- nu0
    for (k in 1:T) {
      v <- v + q * exp(link$log_Delta[, k]) *
        link$D_over_Delta[, k]^2
    }
    I_bb <- crossprod(X, v * X)  ## == X^T diag(v) X, numerically stable since v >= 0.
  } else {
    I_bb <- matrix(0, 0, 0)
  }

  ## ----- α-β cross-block -----
  if (p_beta > 0L) {
    I_ab <- matrix(0, nrow = T, ncol = p_beta)
    ## t = 1: w_{i,1} = q u_{i,1} D_{i,1}/Δ_{i,1} + ν_{i,0}
    w1 <- q * u[, 1] * link$D_over_Delta[, 1] + nu0
    I_ab[1, ] <- as.numeric(crossprod(X, w1))
    if (T >= 2L) {
      for (t in 2:T) {
        w_t <- q * u[, t] * link$D_over_Delta[, t] -
               q * u[, t] * link$D_over_Delta[, t - 1L]
        I_ab[t, ] <- as.numeric(crossprod(X, w_t))
      }
    }
  } else {
    I_ab <- matrix(0, T, 0)
  }

  ## Assemble.
  m <- T + p_beta
  I <- matrix(0, nrow = m, ncol = m)
  I[1:T, 1:T] <- I_aa
  if (p_beta > 0L) {
    I[1:T, (T + 1):m] <- I_ab
    I[(T + 1):m, 1:T] <- t(I_ab)
    I[(T + 1):m, (T + 1):m] <- I_bb
  }
  I
}


## Information matrix in θ = (ã, β): I_θ = J̃^T I_(α,β) J̃.
infor_mat_cumu <- function(theta, X, q, T) {
  atilde <- theta[1:T]
  alpha <- alpha_from_atilde(atilde)
  beta <- if (length(theta) > T) theta[(T + 1):length(theta)] else numeric(0)
  I_ab <- infor_mat_alpha_beta(alpha, beta, X, q)
  J <- full_jacobian(atilde, length(beta))
  crossprod(J, I_ab) %*% J
}


## --- warm-start from data ------------------------------------------------

## §7 initialisation. Returns θ = (ã, β = 0) of length T + p_beta.
warm_start_theta <- function(M, q, T, p_beta, eps = 1e-3) {
  validate_cumu_response(M, T)
  q_mean <- mean(q)
  alpha_hat <- numeric(T)
  for (t in 1:T) {
    rate <- mean(M >= t) / max(q_mean, eps)
    rate <- min(max(rate, eps), 1 - eps)
    alpha_hat[t] <- log(rate / (1 - rate))
  }
  ## Enforce strict ordering before mapping to ã.
  if (T >= 2L) {
    for (t in 2:T) {
      if (alpha_hat[t] >= alpha_hat[t - 1] - eps) {
        alpha_hat[t] <- alpha_hat[t - 1] - eps
      }
    }
  }
  atilde <- numeric(T)
  atilde[1] <- alpha_hat[1]
  if (T >= 2L) {
    for (t in 2:T) {
      atilde[t] <- log(max(eps, alpha_hat[t - 1] - alpha_hat[t]))
    }
  }
  c(atilde, rep.int(0, p_beta))
}


## --- Fisher scoring (unpenalized MLE on ã, β) ----------------------------

## Mirrors irls_iter() in R/param-estimate_logit_get_p_by_t_June.R: loops
## over scoring steps using I^{-1} score_lik. Step-halving ensures that every
## accepted update has a finite, non-decreasing log-likelihood.
## Convergence semantics: 1 = converged, 2 = singular/non-finite system,
## 3 = max-iter, 4 = no acceptable step.
irls_iter_cumu <- function(M_vec, X, theta_estimated, q_vec, T,
                           tolerance = .Machine$double.eps^0.5,
                           stop_criteria = 1e-6, max_iter = 25L,
                           max_halving = 25L) {
  validate_cumu_response(M_vec, T)
  indi <- 1
  n_iter <- 1L
  conv_stat <- 1L
  current_ll <- loss_fun_cumu(theta_estimated, X, M_vec, q_vec, T)
  if (!is.finite(current_ll)) {
    return(c(as.vector(theta_estimated), 4L))
  }
  while (indi >= stop_criteria && n_iter <= max_iter) {
    inf_mat <- infor_mat_cumu(theta_estimated, X, q_vec, T)
    if (any(!is.finite(inf_mat))) {
      conv_stat <- 2L
      break
    }
    ld <- determinant.matrix(inf_mat, logarithm = TRUE)
    if (ld$sign[1] <= 0 || as.numeric(ld$modulus) < log(tolerance)) {
      conv_stat <- 2L
      break
    }
    s_lik <- loss_gradient_cumu(theta_estimated, X, M_vec, q_vec, T)
    if (any(!is.finite(s_lik))) {
      conv_stat <- 2L
      break
    }
    inf_inv <- try(solve(inf_mat), silent = TRUE)
    if (is_error_cumu(inf_inv) || any(!is.finite(inf_inv))) {
      conv_stat <- 2L
      break
    }
    update <- as.numeric(inf_inv %*% s_lik)
    if (any(!is.finite(update))) {
      conv_stat <- 2L
      break
    }
    step_scale <- 1
    accepted <- FALSE
    for (halving in 0:max_halving) {
      theta_update <- theta_estimated + step_scale * update
      proposed_ll <- loss_fun_cumu(theta_update, X, M_vec, q_vec, T)
      if (is.finite(proposed_ll) && proposed_ll >= current_ll - 1e-10) {
        accepted <- TRUE
        break
      }
      step_scale <- step_scale / 2
    }
    if (!accepted) {
      conv_stat <- 4L
      break
    }
    indi <- sum((theta_update - theta_estimated)^2)
    theta_estimated <- theta_update
    current_ll <- proposed_ll
    n_iter <- n_iter + 1L
  }
  if (conv_stat == 1L && indi >= stop_criteria) conv_stat <- 3L
  c(as.vector(theta_estimated), conv_stat)
}


## Null-fit version: hold a subset of θ at zero (mirrors irls_iter_null).
irls_iter_cumu_null <- function(M_vec, X, theta_estimated, hold_zero,
                                q_vec, T,
                                tolerance = .Machine$double.eps^0.5,
                                stop_criteria = 1e-6, max_iter = 25L,
                                max_halving = 25L) {
  validate_cumu_response(M_vec, T)
  indi <- 1
  n_iter <- 1L
  conv_stat <- 1L
  poi <- setdiff(seq_along(theta_estimated), hold_zero)
  current_ll <- loss_fun_cumu(theta_estimated, X, M_vec, q_vec, T)
  if (!is.finite(current_ll)) {
    return(c(as.vector(theta_estimated), 4L))
  }
  while (indi >= stop_criteria && n_iter <= max_iter) {
    inf_mat <- infor_mat_cumu(theta_estimated, X, q_vec, T)
    if (any(!is.finite(inf_mat))) {
      conv_stat <- 2L
      break
    }
    ld <- determinant.matrix(inf_mat[poi, poi, drop = FALSE], logarithm = TRUE)
    if (ld$sign[1] <= 0 || as.numeric(ld$modulus) < log(tolerance)) {
      conv_stat <- 2L
      break
    }
    s_lik <- loss_gradient_cumu(theta_estimated, X, M_vec, q_vec, T)
    if (any(!is.finite(s_lik))) {
      conv_stat <- 2L
      break
    }
    inf_inv <- try(solve(inf_mat[poi, poi, drop = FALSE]), silent = TRUE)
    if (is_error_cumu(inf_inv) || any(!is.finite(inf_inv))) {
      conv_stat <- 2L
      break
    }
    update_poi <- as.numeric(inf_inv %*% s_lik[poi])
    if (any(!is.finite(update_poi))) {
      conv_stat <- 2L
      break
    }
    update_full <- numeric(length(theta_estimated))
    update_full[poi] <- update_poi
    update_full[hold_zero] <- 0
    step_scale <- 1
    accepted <- FALSE
    for (halving in 0:max_halving) {
      theta_update <- theta_estimated + step_scale * update_full
      theta_update[hold_zero] <- 0
      proposed_ll <- loss_fun_cumu(theta_update, X, M_vec, q_vec, T)
      if (is.finite(proposed_ll) && proposed_ll >= current_ll - 1e-10) {
        accepted <- TRUE
        break
      }
      step_scale <- step_scale / 2
    }
    if (!accepted) {
      conv_stat <- 4L
      break
    }
    indi <- sum((theta_update - theta_estimated)^2)
    theta_estimated <- theta_update
    current_ll <- proposed_ll
    n_iter <- n_iter + 1L
  }
  if (conv_stat == 1L && indi >= stop_criteria) conv_stat <- 3L
  c(as.vector(theta_estimated), conv_stat)
}


## --- per-feature drivers --------------------------------------------------

## Mirrors estimate_parameters() shape: takes a region-by-cell matrix and
## fits irls_iter_cumu per region. Returns (T+p+1) × n_features matrix
## where the last row is the convergence code.
estimate_parameters_cumu <- function(r_by_c, X, theta_initial, q_vec, T,
                                     mc.cores = 1L) {
  if ("sparseMatrix" %in% is(r_by_c)) {
    r_by_c <- as.matrix(r_by_c)
  }
  cells_by_region <- t(r_by_c)
  cells_by_region <- as.list(as.data.frame(cells_by_region))
  ## Closure avoids name collision between mclapply's first parameter `X`
  ## and the design matrix `X` we need to pass to irls_iter_cumu.
  fit_one <- function(M_vec) {
    irls_iter_cumu(M_vec, X = X, theta_estimated = theta_initial,
                   q_vec = q_vec, T = T)
  }
  res <- parallel::mclapply(cells_by_region, FUN = fit_one,
                            mc.cores = mc.cores)
  do.call(cbind, res)
}


estimate_parameters_cumu_null <- function(r_by_c, X, theta_initial, hold_zero,
                                          q_vec, T, mc.cores = 1L) {
  if ("sparseMatrix" %in% is(r_by_c)) {
    r_by_c <- as.matrix(r_by_c)
  }
  cells_by_region <- t(r_by_c)
  cells_by_region <- as.list(as.data.frame(cells_by_region))
  fit_one <- function(M_vec) {
    irls_iter_cumu_null(M_vec, X = X, theta_estimated = theta_initial,
                        hold_zero = hold_zero, q_vec = q_vec, T = T)
  }
  res <- parallel::mclapply(cells_by_region, FUN = fit_one,
                            mc.cores = mc.cores)
  do.call(cbind, res)
}
