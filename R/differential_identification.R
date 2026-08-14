### differential_identification


## Ordinary likelihood-ratio test under the exact cumulative-logit model
## (method = "exact" in pacs_test_cumu). Fitting routines live in
## R/param_estimate_cumu.R and the math is in notes/cumulative_logit_math.md
## (§9). Inference is returned only for converged, interior fits.
##
## Args:
##   x_full              : n × p design matrix (alpha-free; alpha lives in theta).
##   theta_estimated_full: (T + p) × n_features matrix of fitted thetas (full).
##   theta_estimated_null: (T + p) × n_features matrix of fitted thetas (null).
##   q_vec               : length-n capture rates.
##   c_by_r              : n × n_features observed-count matrix.
##   T                   : number of thresholds.
##   df_test             : degrees of freedom = number of beta components held
##                         to zero under null (thresholds are shared).
##   conv_full/conv_null : optional per-feature convergence codes. Code 1 is
##                         the only status eligible for inference.
##   mc.cores            : passed to mclapply.
compare_models_cumu <- function(x_full, theta_estimated_full,
                                theta_estimated_null,
                                q_vec, c_by_r, T, df_test,
                                conv_full = NULL, conv_null = NULL,
                                boundary_eps = 1e-3,
                                mc.cores = 1L) {
  if ("sparseMatrix" %in% is(c_by_r)) {
    c_by_r <- as.matrix(c_by_r)
  }
  n_features <- ncol(theta_estimated_full)
  if (n_features != ncol(c_by_r)) {
    stop("theta dimension does not match c_by_r")
  }
  if (ncol(theta_estimated_null) != n_features) {
    stop("full and null theta matrices must have the same number of features")
  }
  expected_rows <- T + ncol(x_full)
  if (nrow(theta_estimated_full) != expected_rows ||
      nrow(theta_estimated_null) != expected_rows) {
    stop(sprintf(
      "theta matrices must each have T + ncol(x_full) = %d rows",
      expected_rows
    ))
  }
  if (nrow(c_by_r) != nrow(x_full) || length(q_vec) != nrow(x_full)) {
    stop("x_full, c_by_r, and q_vec must describe the same cells")
  }
  if (is.null(conv_full)) conv_full <- rep.int(1L, n_features)
  if (is.null(conv_null)) conv_null <- rep.int(1L, n_features)
  if (length(conv_full) != n_features || length(conv_null) != n_features) {
    stop("convergence-code vectors must have one value per feature")
  }

  finite_theta <- apply(is.finite(theta_estimated_full), 2, all) &
    apply(is.finite(theta_estimated_null), 2, all)
  converged <- !is.na(conv_full) & !is.na(conv_null) &
    conv_full == 1L & conv_null == 1L & finite_theta

  ## A standard chi-square reference is not justified when either fit is on
  ## the order-constraint boundary. Flag both models and withhold the p-value.
  boundary <- rep.int(FALSE, n_features)
  if (T >= 2L) {
    near_boundary <- function(theta_mat) {
      colSums(
        exp(theta_mat[2:T, , drop = FALSE]) < boundary_eps,
        na.rm = TRUE
      ) > 0L
    }
    boundary <- near_boundary(theta_estimated_full) |
      near_boundary(theta_estimated_null)
  }
  boundary_eligible <- boundary & converged

  per_peak <- function(j) {
    if (!converged[j] || boundary_eligible[j]) {
      return(list(p = NA_real_, negative = FALSE))
    }
    M_j <- as.numeric(c_by_r[, j])
    th_full <- theta_estimated_full[, j]
    th_null <- theta_estimated_null[, j]
    ll_full <- try(loss_fun_cumu(th_full, x_full, M_j, q_vec, T), silent = TRUE)
    ll_null <- try(loss_fun_cumu(th_null, x_full, M_j, q_vec, T), silent = TRUE)
    if (inherits(ll_full, "try-error") || inherits(ll_null, "try-error") ||
        !is.finite(ll_full) || !is.finite(ll_null)) {
      return(list(p = NA_real_, negative = FALSE))
    }
    stat <- 2 * (ll_full - ll_null)
    numerical_tol <- 1e-8 * (1 + abs(ll_full) + abs(ll_null))
    if (stat < -numerical_tol) {
      return(list(p = NA_real_, negative = TRUE))
    }
    stat <- max(0, stat)
    list(
      p = pchisq(stat, df = df_test, lower.tail = FALSE),
      negative = FALSE
    )
  }

  results <- parallel::mclapply(seq_len(n_features), per_peak,
                                mc.cores = mc.cores)
  pvals <- vapply(results, `[[`, numeric(1L), "p")
  materially_negative <- vapply(results, `[[`, logical(1L), "negative")

  if (any(!converged)) {
    warning(sprintf(
      "%d peak(s) had a non-converged null or full fit; p-values were set to NA.",
      sum(!converged)
    ), call. = FALSE)
  }
  if (any(boundary_eligible)) {
    warning(sprintf(
      paste0(
        "%d peak(s) had a null or full fit on the order-constraint ",
        "boundary (threshold gap < %g); p-values were set to NA."
      ),
      sum(boundary_eligible), boundary_eps
    ), call. = FALSE)
  }
  if (any(materially_negative)) {
    warning(sprintf(
      "%d peak(s) had a materially negative likelihood-ratio statistic; p-values were set to NA.",
      sum(materially_negative)
    ), call. = FALSE)
  }

  pvals
}


#' Compute loss function with Firth regularization from Wii matrix
#'
#' @param wii_sqrt The square root of Wii diagonal elements as a vector
#' @param xdumm Design matrix
#'
#' @return Regularized log loss value
#' @export
#'
loss_firth_from_wii <- function(wii_sqrt, xdumm) {
  wii_sqrt_X <- wii_sqrt * xdumm
  inf_mat <- crossprod(wii_sqrt_X)
  log_loss_star <- 0.5 * determinant.matrix(inf_mat)$modulus[1]
  return(log_loss_star)
}



#' Use gLRT to compare two models and obtain p values
#'
#' @importFrom tictoc tic
#' @importFrom tictoc toc
#' @importFrom parallel mclapply
#' @importFrom stats pchisq
#' @importFrom methods is
#'
#' @param x_full Design matrix for the full model
#' @param theta_estimated_full Estimated coefficients for the full model
#' @param x_null Design matrix for the null model
#' @param theta_estimated_null Estimated coefficients for the null model
#' @param q_vec A vector of capturing probability for each cell
#' @param c_by_r Cell by region matrix
#' @param df_test Degree of freedom of the test
#' @param mc.cores Number of cores in multi-core computing
#'
#' @return A vector of p values based on gLRT
#' @export
compare_models <- function(x_full, theta_estimated_full,
                           x_null, theta_estimated_null,
                           q_vec, c_by_r, df_test, mc.cores = 1) {
  if ("sparseMatrix" %in% is(c_by_r)) {
    c_by_r <- as.matrix(c_by_r)
  }

  ## get number of regions (peaks)
  n_regions <- dim(theta_estimated_full)[2]

  if (n_regions != dim(c_by_r)[2]) {
    stop("please make sure the input matrix and
         the parameters match dimensions")
  }

  ## df
  if (is.vector(theta_estimated_null)) {
    theta_estimated_null <- matrix(theta_estimated_null, nrow = 1)
  }

  ###############
  ### for full matrix
  ###############

  ## vectorize -- calculate for all regions
  x_times_theta <- x_full %*% theta_estimated_full ## n by r (#features)
  p_bg <- 1 - 1 / (exp(x_times_theta) + 1) ## n by r (#features)

  ## calculate loss without Firth prior
  log_pq <- log(q_vec) + log(p_bg) ## n by r (number of features)
  log_1_pq <- log(1 - p_bg * q_vec) ## n by r (number of features)
  log_loss <- colSums(c_by_r * log_pq + (1 - c_by_r) * log_1_pq) ## vector, r

  ## calculate W for GLM
  wii <- (p_bg * (1 - p_bg) * (1 - p_bg) * q_vec) / (1 - p_bg * q_vec) ## n by r
  wii_sqrt <- sqrt(wii) ## n by r

  ###############
  ### for null matrix
  ###############

  ## vectorize -- calculate for all regions
  x_times_theta_null <- x_null %*% theta_estimated_null ## n by r (# features)
  p_bg_null <- 1 - 1 / (exp(x_times_theta_null) + 1) ## n by r (#features)

  ## calculate loss without Firth prior
  log_pq_null <- log(q_vec) + log(p_bg_null) ## n by r (number of features)
  log_1_pq_null <- log(1 - p_bg_null * q_vec) ## n by r (number of features)
  log_loss_null <- colSums(c_by_r * log_pq_null +
    (1 - c_by_r) * log_1_pq_null) ## vector, r

  ## calculate W for GLM
  wii_null <- (p_bg_null * (1 - p_bg_null) * (1 - p_bg_null) * q_vec) /
    (1 - p_bg_null * q_vec) ## n by r
  wii_sqrt_null <- sqrt(wii_null) ## n by r

  ## for each column, do this
  tictoc::tic("computing firth loss start")
  log_loss_firth_full <- mclapply(
    as.data.frame(wii_sqrt), loss_firth_from_wii,
    xdumm = x_full, mc.cores = mc.cores
  )
  log_loss_firth_null <- mclapply(
    as.data.frame(wii_sqrt_null),
    loss_firth_from_wii,
    xdumm = x_null,
    mc.cores = mc.cores
  )
  tictoc::toc()

  log_loss_star_full <- unlist(log_loss_firth_full) + log_loss ## vector len r
  log_loss_star_null <- unlist(log_loss_firth_null) + log_loss_null ## len r

  test_stat <- 2 * (log_loss_star_full - log_loss_star_null)
  LRT_p_value <- pchisq(test_stat, df = df_test, lower.tail = FALSE)

  return(LRT_p_value)
}





#' Internal evaluations only -- Obtain p values using PACS firth model
#'
#' @param para_est estimated coefficients values
#' @param para_reduced_est estimated coefficients values for reduced model
#' @param c_by_r_full Cell by region count matrix
#' @param x_full Design matrix for the full model
#' @param x_null Design matrix for the alternative model
#' @param q_vec A vector of capturing probability
#' @param df_test degree of freedom of the test
#'
#' @return p values
#' @noRd
get_our_firth_p_value <- function(
    para_est, para_reduced_est,
    c_by_r_full, x_full, x_null, q_vec, df_test) {
  lrt_p_val <- vector(length = dim(c_by_r_full)[2])
  names(lrt_p_val) <- colnames(c_by_r_full)

  n_divide <- ceiling(dim(c_by_r_full)[2] / 5000)

  for (i in 1:n_divide) {
    from_i <- 1 + (i - 1) * 5000
    to_i <- min(i * 5000, dim(c_by_r_full)[2])

    c_by_r <- as.matrix(c_by_r_full[, from_i:to_i])
    features <- colnames(c_by_r)

    lrt_p_val[features] <- compare_models(
      x_full = x_full,
      theta_estimated_full = para_est[, from_i:to_i],
      x_null = x_null,
      theta_estimated_null = para_reduced_est[, from_i:to_i],
      q_vec = q_vec,
      c_by_r = c_by_r,
      df_test = df_test,
      mc.cores = 1
    )
  }
  return(lrt_p_val)
}




#' Internal evaluation use only -- loss function Without Firth
#'
#' @param p_vec A vector of open probability
#' @param q_vec A vector of capturing probability
#' @param y_mat Read count matrix
#'
#' @return log loss values
#' @noRd
loss_fun_simple_pq <- function(p_vec, q_vec, y_mat) {
  if (!is.matrix(y_mat)) {
    y_mat <- as.matrix(y_mat)
  }

  p_q_mat <- outer(p_vec, q_vec, FUN = "*")
  log_p <- log(p_vec)
  log_q <- log(q_vec)
  log_pq <- outer(log_p, log_q, FUN = "+") ## better numerical accuracy
  log_1_pq <- log(1 - p_q_mat)

  return(rowSums(y_mat * log_pq + (1 - y_mat) * log_1_pq))
}



#' Internal evaluation use only -- PACS Without Firth
#'
#' @importFrom methods is
#' @param data_matrix_pos Count matrix in foreground group
#' @param data_matrix_neg Count matrix in background group
#' @param true_q_pos Capturing probability vector for the foreground group
#' @param true_q_neg Capturing probability vector for the background group
#'
#' @return P values for each cell
#' @export
our_method_no_firth <- function(
    data_matrix_pos, data_matrix_neg, true_q_pos, true_q_neg) {
  if (!is.matrix(data_matrix_pos)) {
    data_matrix_pos <- as.matrix(data_matrix_pos)
    data_matrix_neg <- as.matrix(data_matrix_neg)
  }

  both_group <- cbind(data_matrix_pos, data_matrix_neg)
  q_truth <- c(true_q_pos, true_q_neg)
  p0_est <- rowSums(both_group) / sum(q_truth)
  p0_est[p0_est > 0.999] <- 0.999 ## make sure it does not exceed 1
  p0_est[p0_est < 0.0001] <- 0.0001 ## make sure it does not exceed 1

  p0_est_A <- rowSums(data_matrix_pos) / sum(true_q_pos)
  p0_est_A[p0_est_A > 0.999] <- 0.999 ## make sure it does not exceed 1
  p0_est_A[p0_est_A < 0.0001] <- 0.0001 ## make sure it does not exceed 1
  p0_est_B <- rowSums(data_matrix_neg) / sum(true_q_neg)
  p0_est_B[p0_est_B > 0.999] <- 0.999 ## make sure it does not exceed 1
  p0_est_B[p0_est_B < 0.0001] <- 0.0001 ## make sure it does not exceed 1

  loss_full_model2 <- loss_fun_simple_pq(
    p_vec = p0_est_A,
    q_vec = true_q_pos,
    y_mat = data_matrix_pos
  ) +
    loss_fun_simple_pq(
      p_vec = p0_est_B,
      q_vec = true_q_neg,
      y_mat = data_matrix_neg
    )

  loss_null_model2 <- loss_fun_simple_pq(
    p_vec = p0_est,
    q_vec = q_truth,
    y_mat = both_group
  )

  LR_value2 <- 2 * (loss_full_model2 - loss_null_model2)
  LR_value2[LR_value2 < 0] <- 0
  sss <- pchisq(LR_value2, df = 1, lower.tail = FALSE)
  return(sss)
}
