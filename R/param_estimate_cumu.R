## Cumulative-logit PACS, exact (Option B)
##
## Implements the proper proportional-odds likelihood with cell-level
## all-or-nothing capture (notes/cumulative_logit_math.md, §3.1, §4.1, §5.1,
## §6, §7). Internal helpers; the public surface is pacs_test_cumu.
##
## At T = 1 these functions reduce term-for-term to loss_fun, loss_gradient,
## infor_mat in R/param-estimate_logit_get_p_by_t_June.R; tests in
## tests/testthat/test-cumu-T1-parity.R verify this.


is_error_cumu <- function(x) inherits(x, "try-error")


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
## Returns a list with n×T matrices p, u, Delta, D.
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
  p <- 1 / (1 + exp(-eta))
  u <- p * (1 - p)
  ## Delta[i, t] = p[i, t] - p[i, t+1] with p[, T+1] = 0.
  if (T >= 2L) {
    Delta <- cbind(p[, -T, drop = FALSE] - p[, -1, drop = FALSE], p[, T])
  } else {
    Delta <- p
  }
  ## D[i, t] = u[i, t] - u[i, t+1] with u[, T+1] = 0.
  if (T >= 2L) {
    D <- cbind(u[, -T, drop = FALSE] - u[, -1, drop = FALSE], u[, T])
  } else {
    D <- u
  }
  list(p = p, u = u, Delta = Delta, D = D)
}


## --- log-likelihood (Option B) -------------------------------------------

## Per-cell log Pr(M_i | x_i, q_i) under Option B (note §3.1). Drops the
## constant log q_i term for m_i >= 1 (irrelevant to optimisation).
loss_fun_cumu <- function(theta, X, M, q, T) {
  alpha <- alpha_from_atilde(theta[1:T])
  beta <- if (length(theta) > T) theta[(T + 1):length(theta)] else numeric(0)
  L <- cumu_link(alpha, beta, X)
  Delta <- L$Delta
  p1 <- L$p[, 1]
  ll <- numeric(length(M))
  zero <- (M == 0L)
  if (any(zero)) {
    ll[zero] <- log(1 - q[zero] * p1[zero])
  }
  if (any(!zero)) {
    idx <- which(!zero)
    delta_picked <- Delta[cbind(idx, M[idx])]
    ll[idx] <- log(delta_picked)
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
  n <- length(M)
  if (is.null(link)) link <- cumu_link(alpha, beta, X)
  p1 <- link$p[, 1]
  u  <- link$u
  Delta <- link$Delta
  D <- link$D

  s_alpha <- numeric(T)
  s_beta <- numeric(p_beta)

  zero <- (M == 0L)
  ## ----- m_i = 0 -----
  if (any(zero)) {
    iz <- which(zero)
    coef_zero <- -q[iz] * u[iz, 1] / (1 - q[iz] * p1[iz])
    s_alpha[1] <- s_alpha[1] + sum(coef_zero)
    if (p_beta > 0L) {
      s_beta <- s_beta + as.numeric(crossprod(X[iz, , drop = FALSE], coef_zero))
    }
  }
  ## ----- m_i >= 1 -----
  if (any(!zero)) {
    inz <- which(!zero)
    Mi <- M[inz]
    delta_pick <- Delta[cbind(inz, Mi)]
    ## ∂ℓ/∂β = (D_{i,M_i} / Delta_{i,M_i}) x_i
    D_pick <- D[cbind(inz, Mi)]
    coef_beta <- D_pick / delta_pick
    if (p_beta > 0L) {
      s_beta <- s_beta + as.numeric(crossprod(X[inz, , drop = FALSE], coef_beta))
    }
    ## ∂ℓ/∂α_t = u_{it}/Delta_{i,M_i} · 1[t=M_i] - u_{it}/Delta_{i,M_i} · 1[t=M_i+1]
    ## Group cells by Mi for vectorised accumulation.
    for (m in unique(Mi)) {
      sub <- which(Mi == m)
      cells <- inz[sub]
      delt <- delta_pick[sub]
      ## t = m contribution (positive)
      s_alpha[m] <- s_alpha[m] + sum(u[cells, m] / delt)
      ## t = m + 1 contribution (negative), only if m + 1 <= T
      if ((m + 1L) <= T) {
        s_alpha[m + 1L] <- s_alpha[m + 1L] - sum(u[cells, m + 1L] / delt)
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
  n <- nrow(X)
  if (is.null(link)) link <- cumu_link(alpha, beta, X)
  p <- link$p; u <- link$u; Delta <- link$Delta; D <- link$D
  p1 <- p[, 1]

  ## ν_{i,0} = q^2 u_{i,1}^2 / (1 - q p_{i,1}); contribution from m=0.
  nu0 <- q^2 * u[, 1]^2 / (1 - q * p1)

  ## ----- α-block (tridiagonal) -----
  I_aa <- matrix(0, nrow = T, ncol = T)
  ## diagonal
  I_aa[1, 1] <- sum(q * u[, 1]^2 / Delta[, 1]) + sum(nu0)
  if (T >= 2L) {
    for (t in 2:T) {
      I_aa[t, t] <- sum(q * u[, t]^2 / Delta[, t]) +
                    sum(q * u[, t]^2 / Delta[, t - 1])
    }
    ## off-diagonal (super-diagonal; symmetric copy below)
    for (t in 1:(T - 1L)) {
      val <- -sum(q * u[, t] * u[, t + 1] / Delta[, t])
      I_aa[t, t + 1L] <- val
      I_aa[t + 1L, t] <- val
    }
  }

  ## ----- β-block -----
  if (p_beta > 0L) {
    v <- nu0
    for (k in 1:T) {
      v <- v + q * D[, k]^2 / Delta[, k]
    }
    I_bb <- crossprod(X, v * X)  ## == X^T diag(v) X, numerically stable since v >= 0.
  } else {
    I_bb <- matrix(0, 0, 0)
  }

  ## ----- α-β cross-block -----
  if (p_beta > 0L) {
    I_ab <- matrix(0, nrow = T, ncol = p_beta)
    ## t = 1: w_{i,1} = q u_{i,1} D_{i,1}/Δ_{i,1} + ν_{i,0}
    w1 <- q * u[, 1] * D[, 1] / Delta[, 1] + nu0
    I_ab[1, ] <- as.numeric(crossprod(X, w1))
    if (T >= 2L) {
      for (t in 2:T) {
        w_t <- q * u[, t] * D[, t] / Delta[, t] -
               q * u[, t] * D[, t - 1] / Delta[, t - 1]
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


## --- Firth penalty (numerical FD on ∂I/∂θ_r) -----------------------------

## Penalty score: 0.5 * tr(I^{-1} ∂I/∂θ_r), r = 1..(T + p_beta).
##
## Phase 1 uses central finite differences on infor_mat_cumu, which costs
## 2*(T + p_beta) Fisher-info evaluations per IRLS step. This is the
## dominant cost of the inner loop and is responsible for the loosened
## 1e-4 tolerance in the T=1 IRLS parity test (vs. 1e-10 for the
## individual loss/score/info parity). An analytic ∂I/∂θ_r is in
## scope for a follow-up using ∂u_{it}/∂η = u_{it}(1 − 2 p_{it}); see
## the math note §5 and follow-up issue tracking the replacement.
loss_grad_pen_cumu <- function(theta, X, q, T,
                               inf_mat = NULL, h = 1e-5) {
  if (is.null(inf_mat)) inf_mat <- infor_mat_cumu(theta, X, q, T)
  m <- length(theta)
  i_inv <- try(solve(inf_mat), silent = TRUE)
  if (is_error_cumu(i_inv)) return(rep.int(NA_real_, m))
  out <- numeric(m)
  for (r in seq_len(m)) {
    th_p <- theta; th_p[r] <- th_p[r] + h
    th_m <- theta; th_m[r] <- th_m[r] - h
    I_p <- infor_mat_cumu(th_p, X, q, T)
    I_m <- infor_mat_cumu(th_m, X, q, T)
    dI <- (I_p - I_m) / (2 * h)
    out[r] <- 0.5 * sum(diag(i_inv %*% dI))
  }
  out
}


## Penalised log-likelihood loss + 0.5 log det I.
loss_fun_star_cumu <- function(theta, X, M, q, T, inf_mat = NULL) {
  ll <- loss_fun_cumu(theta, X, M, q, T)
  if (is.null(inf_mat)) inf_mat <- infor_mat_cumu(theta, X, q, T)
  ld <- determinant.matrix(inf_mat, logarithm = TRUE)
  if (ld$sign[1] <= 0) return(-Inf)
  ll + 0.5 * as.numeric(ld$modulus)
}


## --- Barndorff-Nielsen r* saddlepoint adjustment (df = 1) ----------------

## EXPERIMENTAL. For testing a scalar ψ (one β component), replaces the
## first-order chi-squared tail probability with the Barndorff-Nielsen r*
## approximation, which targets O(n^{-3/2}) tail accuracy.
##
## Important: simulation shows the Firth-corrected chi-squared reference
## is already very well calibrated (KS > 0.4 at n=300, > 0.9 at n=80),
## and the r* correction can DEGRADE calibration because Firth's penalty
## already absorbs the O(1/n) bias that r* targets. The two corrections
## are alternative approaches to the same asymptotic deficiency; combining
## them produces over-correction. Use this only for research comparisons,
## not as a production default. See math note §9.1 for the full analysis.
##
## Implementation uses the unpenalized signed-root LRT with the full
## profile score (no S_λ=0 shortcut, since the Firth MLE zeros ∂ℓ*/∂λ
## not ∂ℓ/∂λ). Falls back to pchisq when: df > 1, r ≈ 0, sign mismatch,
## or non-positive profile information.
saddlepoint_pvalue_scalar <- function(stat_pen, th_full, th_null, psi_idx,
                                      X, M, q, T) {
  if (stat_pen <= 0) return(1.0)

  lam_idx <- setdiff(seq_along(th_full), psi_idx)

  ## Unpenalized LRT signed root (consistent with unpenalized score/info).
  ll_full_unpen <- loss_fun_cumu(th_full, X, M, q, T)
  ll_null_unpen <- loss_fun_cumu(th_null, X, M, q, T)
  stat_unpen <- max(0, 2 * (ll_full_unpen - ll_null_unpen))
  r <- sign(th_full[psi_idx] - th_null[psi_idx]) * sqrt(stat_unpen)

  ## Profile Fisher information at full MLE (Schur complement).
  I_full <- infor_mat_cumu(th_full, X, q, T)
  I_pp <- I_full[psi_idx, psi_idx]
  I_pl <- I_full[psi_idx, lam_idx, drop = FALSE]
  I_ll_full <- I_full[lam_idx, lam_idx, drop = FALSE]
  I_ll_full_inv <- try(solve(I_ll_full), silent = TRUE)
  if (is_error_cumu(I_ll_full_inv)) {
    return(pchisq(stat_pen, df = 1L, lower.tail = FALSE))
  }
  J_profile <- as.numeric(I_pp - I_pl %*% I_ll_full_inv %*% t(I_pl))

  if (J_profile <= 0) {
    return(pchisq(stat_pen, df = 1L, lower.tail = FALSE))
  }

  ## Full profile score at null MLE: S_{ψ.λ} = S_ψ - I_ψλ I_λλ^{-1} S_λ.
  ## (No shortcut: Firth MLE zeros ∂ℓ*/∂λ, not ∂ℓ/∂λ.)
  s_unpen <- loss_gradient_cumu(th_null, X, M, q, T)
  if (anyNA(s_unpen)) {
    return(pchisq(stat_pen, df = 1L, lower.tail = FALSE))
  }
  I_null <- infor_mat_cumu(th_null, X, q, T)
  I_ll_null <- I_null[lam_idx, lam_idx, drop = FALSE]
  I_pl_null <- I_null[psi_idx, lam_idx, drop = FALSE]
  I_ll_null_inv <- try(solve(I_ll_null), silent = TRUE)
  if (is_error_cumu(I_ll_null_inv)) {
    return(pchisq(stat_pen, df = 1L, lower.tail = FALSE))
  }
  S_psi_profile <- s_unpen[psi_idx] -
    as.numeric(I_pl_null %*% I_ll_null_inv %*% s_unpen[lam_idx])
  u <- S_psi_profile / sqrt(J_profile)

  if (abs(r) < 1e-7 || (u / r) <= 0) {
    return(pchisq(stat_pen, df = 1L, lower.tail = FALSE))
  }

  r_star <- r + (1 / r) * log(u / r)
  2 * pnorm(-abs(r_star))
}


## --- warm-start from data ------------------------------------------------

## §7 initialisation. Returns θ = (ã, β = 0) of length T + p_beta.
warm_start_theta <- function(M, q, T, p_beta, eps = 1e-3) {
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


## --- IRLS (Fisher scoring on ã, β) ---------------------------------------

## Mirrors irls_iter() in R/param-estimate_logit_get_p_by_t_June.R: loops
## over Newton steps using I^{-1} (score_lik + score_pen), with the same
## convergence semantics (1 = converged, 2 = singular, 3 = max-iter).
irls_iter_cumu <- function(M_vec, X, theta_estimated, q_vec, T,
                           tolerance = .Machine$double.eps^0.5,
                           stop_criteria = 1e-6, max_iter = 25L) {
  indi <- 1
  n_iter <- 1L
  conv_stat <- 1L
  while (indi >= stop_criteria && n_iter <= max_iter) {
    inf_mat <- infor_mat_cumu(theta_estimated, X, q_vec, T)
    ld <- determinant.matrix(inf_mat, logarithm = TRUE)
    if (ld$sign[1] <= 0 || as.numeric(ld$modulus) < log(tolerance)) {
      conv_stat <- 2L
      break
    }
    s_lik <- loss_gradient_cumu(theta_estimated, X, M_vec, q_vec, T)
    s_pen <- loss_grad_pen_cumu(theta_estimated, X, q_vec, T, inf_mat = inf_mat)
    if (anyNA(s_lik) || anyNA(s_pen)) {
      conv_stat <- 2L
      break
    }
    inf_inv <- try(solve(inf_mat), silent = TRUE)
    if (is_error_cumu(inf_inv)) {
      conv_stat <- 2L
      break
    }
    update <- as.numeric(inf_inv %*% (s_lik + s_pen))
    theta_update <- theta_estimated + update
    indi <- sum((theta_update - theta_estimated)^2)
    theta_estimated <- theta_update
    n_iter <- n_iter + 1L
  }
  if (indi >= stop_criteria) conv_stat <- 3L
  c(as.vector(theta_estimated), conv_stat)
}


## Null-fit version: hold a subset of θ at zero (mirrors irls_iter_null).
irls_iter_cumu_null <- function(M_vec, X, theta_estimated, hold_zero,
                                q_vec, T,
                                tolerance = .Machine$double.eps^0.5,
                                stop_criteria = 1e-6, max_iter = 25L) {
  indi <- 1
  n_iter <- 1L
  conv_stat <- 1L
  poi <- setdiff(seq_along(theta_estimated), hold_zero)
  zero_augment <- rep.int(0, length(hold_zero))
  while (indi >= stop_criteria && n_iter <= max_iter) {
    inf_mat <- infor_mat_cumu(theta_estimated, X, q_vec, T)
    ld <- determinant.matrix(inf_mat[poi, poi, drop = FALSE], logarithm = TRUE)
    if (ld$sign[1] <= 0 || as.numeric(ld$modulus) < log(tolerance)) {
      conv_stat <- 2L
      break
    }
    s_lik <- loss_gradient_cumu(theta_estimated, X, M_vec, q_vec, T)
    s_pen <- loss_grad_pen_cumu(theta_estimated, X, q_vec, T, inf_mat = inf_mat)
    if (anyNA(s_lik) || anyNA(s_pen)) {
      conv_stat <- 2L
      break
    }
    s_total <- s_lik + s_pen
    inf_inv <- try(solve(inf_mat[poi, poi, drop = FALSE]), silent = TRUE)
    if (is_error_cumu(inf_inv)) {
      conv_stat <- 2L
      break
    }
    update_poi <- as.numeric(inf_inv %*% s_total[poi])
    update_full <- numeric(length(theta_estimated))
    update_full[poi] <- update_poi
    update_full[hold_zero] <- 0
    theta_update <- theta_estimated + update_full
    indi <- sum((theta_update - theta_estimated)^2)
    theta_estimated <- theta_update
    n_iter <- n_iter + 1L
  }
  if (indi >= stop_criteria) conv_stat <- 3L
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
