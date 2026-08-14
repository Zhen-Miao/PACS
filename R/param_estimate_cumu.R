## Cumulative-logit PACS, exact likelihood
##
## Implements the proper proportional-odds likelihood with capture-rate
## correction (notes/cumulative_logit_math.md, §3-§7). Two capture models are
## supported, selected by the `capture` argument threaded through every
## routine here:
##
##   capture = "B"  cell-level all-or-nothing dropout (note §3.1, §4.1, §5.1)
##   capture = "A"  fragment-level binomial thinning, M_i | Y_i = k ~
##                  Binom(k, q_i) (note §3.2, §4.2, §5.2)
##
## Option A is the natural generalisation of the binary capture model; Option
## B is cheaper and remains the default for backward compatibility. Both are
## unpenalized maximum likelihood: no bias reduction is claimed (note §6).
##
## At T = 1 the two options coincide, and the likelihood, score, and expected
## information reduce term-for-term to loss_fun, loss_gradient, and infor_mat
## in R/param-estimate_logit_get_p_by_t_June.R; tests in
## tests/testthat/test-cumu-T1-parity.R and
## tests/testthat/test-cumu-optionA.R verify this.


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
## counterparts, plus the n×(T+1) category-indexed matrices log_pi and r used
## by Option A. Delta is evaluated without subtracting two rounded logistic
## probabilities.
#' @importFrom stats plogis
#' @noRd
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

  ## Category-indexed views, columns k = 0, 1, ..., T. Option B only ever
  ## touches the single category k = m_i, but Option A mixes over all
  ## k >= m_i, so it needs the k = 0 category as well.
  ##   pi_{ik} = p_{ik} - p_{i,k+1},   pi_{i0} = 1 - p_{i1}
  ##   r_{ik}  = (u_{ik} - u_{i,k+1}) / pi_{ik} = 1 - p_{ik} - p_{i,k+1},
  ## the latter valid at k = 0 (giving -p_{i1}) and k = T (giving 1 - p_{iT})
  ## under the conventions above.
  log_pi <- cbind(log_one_minus_p[, 1L], log_Delta)
  r <- cbind(-p[, 1L], D_over_Delta)

  list(
    p = p, one_minus_p = one_minus_p,
    log_p = log_p, log_one_minus_p = log_one_minus_p,
    u = u, log_u = log_u,
    Delta = Delta, log_Delta = log_Delta,
    D = D, D_over_Delta = D_over_Delta,
    log_pi = log_pi, r = r
  )
}


## --- Option A: fragment-level binomial thinning ---------------------------
##
## M_i | Y_i = k ~ Binomial(k, q_i), so (note §2, §3.2, §4.2, §5.2)
##
##   Pr(M_i = m) = sum_{k >= m} w_{ik}^{(m)} pi_{ik},
##   w_{ik}^{(m)} = C(k, m) q_i^m (1 - q_i)^{k - m}.
##
## Everything below works from log w and log pi so that the mixture stays
## stable when individual category probabilities underflow.


## log w_{ik}^{(m)} as an n × (T+1) × (T+1) array indexed [cell, m+1, k+1].
## Entries with k < m are -Inf. Depends only on (q, T), never on θ, so a fit
## computes it once and reuses it across scoring iterations.
thinning_log_weights <- function(q, T) {
  n <- length(q)
  log_q <- log(q)
  log_one_minus_q <- log1p(-q)
  log_w <- array(-Inf, dim = c(n, T + 1L, T + 1L))
  for (m in 0:T) {
    ## m = 0 must not touch log_q, and k = m must not touch log(1 - q_i):
    ## at q_i = 1 the latter is -Inf and 0 * -Inf is NaN in R.
    captured <- if (m > 0L) m * log_q else numeric(n)
    for (k in m:T) {
      missed <- if (k > m) (k - m) * log_one_minus_q else numeric(n)
      log_w[, m + 1L, k + 1L] <- lchoose(k, m) + captured + missed
    }
  }
  log_w
}


## Row-wise log-sum-exp. Rows whose entries are all -Inf return -Inf.
row_logsumexp <- function(z) {
  row_max <- z[, 1L]
  if (ncol(z) > 1L) {
    for (j in 2:ncol(z)) row_max <- pmax(row_max, z[, j])
  }
  out <- rep.int(-Inf, length(row_max))
  usable <- is.finite(row_max)
  if (any(usable)) {
    shifted <- exp(z[usable, , drop = FALSE] - row_max[usable])
    out[usable] <- row_max[usable] + log(rowSums(shifted))
  }
  out
}


## n × (T+1) matrix of log Pr(M_i = m | x_i, q_i), columns m = 0, ..., T.
optA_log_prob <- function(link, log_w, T) {
  n <- nrow(link$log_pi)
  out <- matrix(-Inf, nrow = n, ncol = T + 1L)
  for (m in 0:T) {
    k_seq <- m:T
    terms <- matrix(log_w[, m + 1L, k_seq + 1L],
                    nrow = n, ncol = length(k_seq)) +
      link$log_pi[, k_seq + 1L, drop = FALSE]
    out[, m + 1L] <- row_logsumexp(terms)
  }
  out
}


## Per-cell Option A score evaluated at the fixed observed category m, for the
## cells in `rows`. Returns the α part as a length(rows) × T matrix and the β
## part as the scalar multiplier c_i in ∂ℓ_i/∂β = c_i x_i.
##
## The α part uses the cancellation
##   γ_{ik} · u_{it} / pi_{ik} = exp(log w_{ik}^{(m)} + log u_{it} - ℓ_i),
## which follows from log γ_{ik} = log w_{ik}^{(m)} + log pi_{ik} - ℓ_i and
## drops pi_{ik} entirely. The ratio u_{it} / pi_{ik} is never formed: it
## overflows when a category probability underflows, even though the product
## with γ_{ik} stays O(1).
##
## Callers must pass only rows with finite ℓ_i; γ is undefined otherwise.
optA_score_pieces <- function(link, log_w, log_prob_m, m, T, rows) {
  n_rows <- length(rows)
  alpha_part <- matrix(0, nrow = n_rows, ncol = T)
  if (n_rows == 0L) {
    return(list(alpha = alpha_part, beta = numeric(0)))
  }
  log_prob_rows <- log_prob_m[rows]
  k_seq <- m:T

  ## γ_{ik} = Pr(Y_i = k | M_i = m), k = m, ..., T.
  gamma <- exp(
    matrix(log_w[rows, m + 1L, k_seq + 1L],
           nrow = n_rows, ncol = length(k_seq)) +
      link$log_pi[rows, k_seq + 1L, drop = FALSE] - log_prob_rows
  )
  ## ∂ℓ_i/∂β = (sum_k γ_{ik} r_{ik}) x_i; |r_{ik}| < 1 so this needs no
  ## log-space treatment.
  beta_part <- rowSums(gamma * link$r[rows, k_seq + 1L, drop = FALSE])

  ## ∂ℓ_i/∂α_t = γ_{it} u_{it}/pi_{it} - γ_{i,t-1} u_{it}/pi_{i,t-1}, where a
  ## term is present only when its category index is at least m (γ_{ik} = 0
  ## for k < m).
  for (t in seq_len(T)) {
    contribution <- numeric(n_rows)
    if (t >= m) {
      contribution <- contribution +
        exp(log_w[rows, m + 1L, t + 1L] + link$log_u[rows, t] -
              log_prob_rows)
    }
    if (t - 1L >= m) {
      contribution <- contribution -
        exp(log_w[rows, m + 1L, t] + link$log_u[rows, t] - log_prob_rows)
    }
    alpha_part[, t] <- contribution
  }

  list(alpha = alpha_part, beta = beta_part)
}


## Option A score in (α, β), summed over cells. Mirrors score_alpha_beta's
## return shape.
score_alpha_beta_optA <- function(link, X, M, q, T, p_beta, log_w = NULL) {
  if (is.null(log_w)) log_w <- thinning_log_weights(q, T)
  log_prob <- optA_log_prob(link, log_w, T)
  s_alpha <- numeric(T)
  s_beta <- numeric(p_beta)
  for (m in sort(unique(M))) {
    cells <- which(M == m)
    pieces <- optA_score_pieces(link, log_w, log_prob[, m + 1L], m, T,
                                rows = cells)
    s_alpha <- s_alpha + colSums(pieces$alpha)
    if (p_beta > 0L) {
      s_beta <- s_beta +
        as.numeric(crossprod(X[cells, , drop = FALSE], pieces$beta))
    }
  }
  c(s_alpha, s_beta)
}


## Expected Fisher information under Option A (note §5.2). There is no
## closed-form block simplification as under Option B — the α-block is not
## tridiagonal, because a cell with M_i = m puts posterior mass on every
## latent k >= m — so the categories are enumerated directly:
##
##   I(θ) = sum_i sum_{m=0}^{T} Pr(M_i = m) s_{i,m} s_{i,m}^T.
infor_mat_alpha_beta_optA <- function(link, X, q, T, p_beta, log_w = NULL) {
  if (is.null(log_w)) log_w <- thinning_log_weights(q, T)
  log_prob <- optA_log_prob(link, log_w, T)
  n_par <- T + p_beta
  I_aa <- matrix(0, nrow = T, ncol = T)
  I_ab <- matrix(0, nrow = T, ncol = p_beta)
  v_beta <- numeric(nrow(link$log_pi))

  for (m in 0:T) {
    log_prob_m <- log_prob[, m + 1L]
    ## Categories that are numerically impossible contribute nothing and
    ## would give an undefined posterior, so they are skipped rather than
    ## multiplied by a zero weight.
    weight_all <- exp(log_prob_m)
    rows <- which(is.finite(log_prob_m) & weight_all > 0)
    if (length(rows) == 0L) next
    weight <- weight_all[rows]
    pieces <- optA_score_pieces(link, log_w, log_prob_m, m, T, rows = rows)
    I_aa <- I_aa + crossprod(pieces$alpha, weight * pieces$alpha)
    if (p_beta > 0L) {
      I_ab <- I_ab + crossprod(
        pieces$alpha, (weight * pieces$beta) * X[rows, , drop = FALSE]
      )
      v_beta[rows] <- v_beta[rows] + weight * pieces$beta^2
    }
  }

  I <- matrix(0, nrow = n_par, ncol = n_par)
  I[1:T, 1:T] <- I_aa
  if (p_beta > 0L) {
    I[1:T, (T + 1):n_par] <- I_ab
    I[(T + 1):n_par, 1:T] <- t(I_ab)
    I[(T + 1):n_par, (T + 1):n_par] <- crossprod(X, v_beta * X)
  }
  I
}


## --- log-likelihood ------------------------------------------------------

## Summed log Pr(M_i | x_i, q_i) over cells.
##
## capture = "B" (note §3.1) drops the constant log q_i term for m_i >= 1,
## which is irrelevant to optimisation. capture = "A" (note §3.2) returns the
## full log-likelihood; its capture terms depend on m_i through the binomial
## weights and cannot be factored out. Both conventions are internally
## consistent, so likelihood-ratio statistics are unaffected — but the two
## values are not comparable to each other.
loss_fun_cumu <- function(theta, X, M, q, T, validate = TRUE,
                          capture = c("B", "A"), log_w = NULL) {
  capture <- match.arg(capture)
  if (validate) validate_cumu_response(M, T)
  alpha <- alpha_from_atilde(theta[1:T])
  beta <- if (length(theta) > T) theta[(T + 1):length(theta)] else numeric(0)
  L <- cumu_link(alpha, beta, X)
  if (capture == "A") {
    if (is.null(log_w)) log_w <- thinning_log_weights(q, T)
    log_prob <- optA_log_prob(L, log_w, T)
    return(sum(log_prob[cbind(seq_along(M), M + 1L)]))
  }
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
score_alpha_beta <- function(alpha, beta, X, M, q, link = NULL,
                             validate = TRUE, capture = c("B", "A"),
                             log_w = NULL) {
  capture <- match.arg(capture)
  T <- length(alpha)
  p_beta <- length(beta)
  if (is.null(link)) link <- cumu_link(alpha, beta, X)
  if (validate) validate_cumu_response(M, T)
  if (capture == "A") {
    return(score_alpha_beta_optA(link, X, M, q, T, p_beta, log_w))
  }
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
loss_gradient_cumu <- function(theta, X, M, q, T, validate = TRUE,
                               capture = c("B", "A"), log_w = NULL) {
  capture <- match.arg(capture)
  atilde <- theta[1:T]
  alpha <- alpha_from_atilde(atilde)
  beta <- if (length(theta) > T) theta[(T + 1):length(theta)] else numeric(0)
  s_ab <- score_alpha_beta(alpha, beta, X, M, q, validate = validate,
                           capture = capture, log_w = log_w)
  J <- full_jacobian(atilde, length(beta))
  as.numeric(crossprod(J, s_ab))
}


## --- Fisher information in (α, β) ----------------------------------------
##
## Option B uses the closed-form blocks of §5.1; Option A enumerates observed
## categories per §5.2.
infor_mat_alpha_beta <- function(alpha, beta, X, q, link = NULL,
                                 capture = c("B", "A"), log_w = NULL) {
  capture <- match.arg(capture)
  T <- length(alpha)
  p_beta <- length(beta)
  if (is.null(link)) link <- cumu_link(alpha, beta, X)
  if (capture == "A") {
    return(infor_mat_alpha_beta_optA(link, X, q, T, p_beta, log_w))
  }
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
infor_mat_cumu <- function(theta, X, q, T, capture = c("B", "A"),
                           log_w = NULL) {
  capture <- match.arg(capture)
  atilde <- theta[1:T]
  alpha <- alpha_from_atilde(atilde)
  beta <- if (length(theta) > T) theta[(T + 1):length(theta)] else numeric(0)
  I_ab <- infor_mat_alpha_beta(alpha, beta, X, q, capture = capture,
                               log_w = log_w)
  J <- full_jacobian(atilde, length(beta))
  crossprod(J, I_ab) %*% J
}


## --- warm-start from data ------------------------------------------------

## §7 initialisation. Returns θ = (ã, β = 0) of length T + p_beta.
##
## The capture correction mean(M >= t) / mean(q) inverts Option B's marginal
## exactly. Under Option A the relation between Pr(M >= t) and Pr(Y >= t) also
## involves the intermediate categories, so the same expression is only a
## rough inverse there. It is used for both because it is a starting value:
## Fisher scoring with step-halving does the rest.
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

## Shared engine for full and constrained-null fits. Step-halving ensures that
## every accepted update has a finite, non-decreasing log-likelihood. A fit is
## converged only when both the accepted step and the free-parameter score are
## small, so a heavily halved but non-stationary step cannot report success.
## Status: 1 = converged, 2 = singular/non-finite scoring system,
## 3 = max-iter, 4 = no acceptable step, 5 = non-finite starting likelihood.
irls_fit_cumu <- function(M_vec, X, theta_estimated, q_vec, T,
                          hold_zero = integer(0),
                          tolerance = .Machine$double.eps^0.5,
                          stop_criteria = 1e-6,
                          score_tolerance = 1e-4,
                          max_iter = 25L, max_halving = 25L,
                          capture = c("B", "A"), log_w = NULL) {
  capture <- match.arg(capture)
  validate_cumu_response(M_vec, T)
  if (length(max_iter) != 1L || max_iter < 1L || max_iter != as.integer(max_iter)) {
    stop("max_iter must be a positive integer.", call. = FALSE)
  }
  if (length(max_halving) != 1L || max_halving < 0L ||
      max_halving != as.integer(max_halving)) {
    stop("max_halving must be a non-negative integer.", call. = FALSE)
  }
  if (!is.finite(stop_criteria) || stop_criteria <= 0 ||
      !is.finite(score_tolerance) || score_tolerance < 0) {
    stop("convergence tolerances must be finite and non-negative.", call. = FALSE)
  }

  hold_zero <- unique(as.integer(hold_zero))
  if (any(hold_zero < 1L | hold_zero > length(theta_estimated))) {
    stop("hold_zero contains an invalid parameter index.", call. = FALSE)
  }
  free <- setdiff(seq_along(theta_estimated), hold_zero)
  if (length(free) == 0L) {
    stop("at least one parameter must remain free.", call. = FALSE)
  }
  theta_estimated[hold_zero] <- 0

  ## The thinning weights depend only on (q, T), so pay for them once instead
  ## of on every likelihood, score, and information evaluation below.
  if (capture == "A" && is.null(log_w)) {
    log_w <- thinning_log_weights(q_vec, T)
  }

  current_ll <- loss_fun_cumu(
    theta_estimated, X, M_vec, q_vec, T, validate = FALSE,
    capture = capture, log_w = log_w
  )
  if (!is.finite(current_ll)) {
    return(c(as.vector(theta_estimated), 5L))
  }

  conv_stat <- 3L
  n_iter <- 1L
  while (n_iter <= max_iter) {
    inf_mat <- infor_mat_cumu(theta_estimated, X, q_vec, T,
                              capture = capture, log_w = log_w)
    if (any(!is.finite(inf_mat))) {
      conv_stat <- 2L
      break
    }
    inf_free <- inf_mat[free, free, drop = FALSE]
    ld <- determinant.matrix(inf_free, logarithm = TRUE)
    if (ld$sign[1] <= 0 || as.numeric(ld$modulus) < log(tolerance)) {
      conv_stat <- 2L
      break
    }

    score <- loss_gradient_cumu(
      theta_estimated, X, M_vec, q_vec, T, validate = FALSE,
      capture = capture, log_w = log_w
    )
    if (any(!is.finite(score))) {
      conv_stat <- 2L
      break
    }
    inf_inv <- try(solve(inf_free), silent = TRUE)
    if (is_error_cumu(inf_inv) || any(!is.finite(inf_inv))) {
      conv_stat <- 2L
      break
    }
    update_free <- as.numeric(inf_inv %*% score[free])
    if (any(!is.finite(update_free))) {
      conv_stat <- 2L
      break
    }
    update <- numeric(length(theta_estimated))
    update[free] <- update_free

    step_scale <- 1
    accepted <- FALSE
    for (halving in seq.int(0L, as.integer(max_halving))) {
      theta_update <- theta_estimated + step_scale * update
      theta_update[hold_zero] <- 0
      proposed_ll <- loss_fun_cumu(
        theta_update, X, M_vec, q_vec, T, validate = FALSE,
        capture = capture, log_w = log_w
      )
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

    step_size <- sum((theta_update - theta_estimated)^2)
    theta_estimated <- theta_update
    current_ll <- proposed_ll
    updated_score <- loss_gradient_cumu(
      theta_estimated, X, M_vec, q_vec, T, validate = FALSE,
      capture = capture, log_w = log_w
    )
    if (any(!is.finite(updated_score))) {
      conv_stat <- 2L
      break
    }
    score_norm <- max(abs(updated_score[free]))
    if (step_size < stop_criteria && score_norm <= score_tolerance) {
      conv_stat <- 1L
      break
    }
    n_iter <- n_iter + 1L
  }

  c(as.vector(theta_estimated), conv_stat)
}


irls_iter_cumu <- function(M_vec, X, theta_estimated, q_vec, T,
                           tolerance = .Machine$double.eps^0.5,
                           stop_criteria = 1e-6,
                           score_tolerance = 1e-4,
                           max_iter = 25L, max_halving = 25L,
                           capture = c("B", "A"), log_w = NULL) {
  irls_fit_cumu(
    M_vec = M_vec, X = X, theta_estimated = theta_estimated,
    q_vec = q_vec, T = T, tolerance = tolerance,
    stop_criteria = stop_criteria, score_tolerance = score_tolerance,
    max_iter = max_iter, max_halving = max_halving,
    capture = match.arg(capture), log_w = log_w
  )
}


## Null-fit version: hold a subset of θ at zero (mirrors irls_iter_null).
irls_iter_cumu_null <- function(M_vec, X, theta_estimated, hold_zero,
                                q_vec, T,
                                tolerance = .Machine$double.eps^0.5,
                                stop_criteria = 1e-6,
                                score_tolerance = 1e-4,
                                max_iter = 25L, max_halving = 25L,
                                capture = c("B", "A"), log_w = NULL) {
  irls_fit_cumu(
    M_vec = M_vec, X = X, theta_estimated = theta_estimated,
    q_vec = q_vec, T = T, hold_zero = hold_zero,
    tolerance = tolerance, stop_criteria = stop_criteria,
    score_tolerance = score_tolerance, max_iter = max_iter,
    max_halving = max_halving, capture = match.arg(capture), log_w = log_w
  )
}


## --- per-feature drivers --------------------------------------------------

## Mirrors estimate_parameters() shape: takes a region-by-cell matrix and
## fits irls_iter_cumu per region. Returns (T+p+1) × n_features matrix
## where the last row is the convergence code.
estimate_parameters_cumu <- function(r_by_c, X, theta_initial, q_vec, T,
                                     mc.cores = 1L, capture = c("B", "A")) {
  capture <- match.arg(capture)
  if ("sparseMatrix" %in% is(r_by_c)) {
    r_by_c <- as.matrix(r_by_c)
  }
  cells_by_region <- t(r_by_c)
  cells_by_region <- as.list(as.data.frame(cells_by_region))
  ## Shared across peaks: the thinning weights depend only on (q, T).
  log_w <- if (capture == "A") thinning_log_weights(q_vec, T) else NULL
  ## Closure avoids name collision between mclapply's first parameter `X`
  ## and the design matrix `X` we need to pass to irls_iter_cumu.
  fit_one <- function(M_vec) {
    irls_iter_cumu(M_vec, X = X, theta_estimated = theta_initial,
                   q_vec = q_vec, T = T, capture = capture, log_w = log_w)
  }
  res <- parallel::mclapply(cells_by_region, FUN = fit_one,
                            mc.cores = mc.cores)
  do.call(cbind, res)
}


estimate_parameters_cumu_null <- function(r_by_c, X, theta_initial, hold_zero,
                                          q_vec, T, mc.cores = 1L,
                                          capture = c("B", "A")) {
  capture <- match.arg(capture)
  if ("sparseMatrix" %in% is(r_by_c)) {
    r_by_c <- as.matrix(r_by_c)
  }
  cells_by_region <- t(r_by_c)
  cells_by_region <- as.list(as.data.frame(cells_by_region))
  log_w <- if (capture == "A") thinning_log_weights(q_vec, T) else NULL
  fit_one <- function(M_vec) {
    irls_iter_cumu_null(M_vec, X = X, theta_estimated = theta_initial,
                        hold_zero = hold_zero, q_vec = q_vec, T = T,
                        capture = capture, log_w = log_w)
  }
  res <- parallel::mclapply(cells_by_region, FUN = fit_one,
                            mc.cores = mc.cores)
  do.call(cbind, res)
}
