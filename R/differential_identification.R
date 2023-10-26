### differential_identification

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
    stop("please make sure the input matrix and the parameters match dimensions")
  }

  ## df
  if (is.vector(theta_estimated_null)) {
    theta_estimated_null <- matrix(theta_estimated_null, nrow = 1)
  }

  ###############
  ### for full matrix
  ###############

  ## vectorize -- calculate for all regions
  x_times_theta <- x_full %*% theta_estimated_full ## n by r (number of features)
  p_bg <- 1 - 1 / (exp(x_times_theta) + 1) ## n by r (number of features)

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
  x_times_theta_null <- x_null %*% theta_estimated_null ## n by r (number of features)
  p_bg_null <- 1 - 1 / (exp(x_times_theta_null) + 1) ## n by r (number of features)

  ## calculate loss without Firth prior
  log_pq_null <- log(q_vec) + log(p_bg_null) ## n by r (number of features)
  log_1_pq_null <- log(1 - p_bg_null * q_vec) ## n by r (number of features)
  log_loss_null <- colSums(c_by_r * log_pq_null + (1 - c_by_r) * log_1_pq_null) ## vector, r

  ## calculate W for GLM
  wii_null <- (p_bg_null * (1 - p_bg_null) * (1 - p_bg_null) * q_vec) / (1 - p_bg_null * q_vec) ## n by r
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

  log_loss_star_full <- unlist(log_loss_firth_full) + log_loss ## vector length r
  log_loss_star_null <- unlist(log_loss_firth_null) + log_loss_null ## vector length r

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
#' @importFrom Rfast rowsums
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
  # print(quantile(p_q_mat))
  log_p <- log(p_vec)
  log_q <- log(q_vec)
  log_pq <- outer(log_p, log_q, FUN = "+") ## better numerical accuracy
  log_1_pq <- log(1 - p_q_mat)

  return(Rfast::rowsums(y_mat * log_pq + (1 - y_mat) * log_1_pq))
}



#' Internal evaluation use only -- PACS Without Firth
#'
#' @importFrom Rfast rowsums
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
  p0_est <- Rfast::rowsums(both_group) / sum(q_truth)
  p0_est[p0_est > 0.999] <- 0.999 ## make sure it does not exceed 1
  p0_est[p0_est < 0.0001] <- 0.0001 ## make sure it does not exceed 1

  p0_est_A <- Rfast::rowsums(data_matrix_pos) / sum(true_q_pos)
  p0_est_A[p0_est_A > 0.999] <- 0.999 ## make sure it does not exceed 1
  p0_est_A[p0_est_A < 0.0001] <- 0.0001 ## make sure it does not exceed 1
  p0_est_B <- Rfast::rowsums(data_matrix_neg) / sum(true_q_neg)
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
  # LR_value2[LR_value2 < 0] = 0
  sss <- pchisq(LR_value2, df = 1, lower.tail = FALSE)
  return(sss)
}
