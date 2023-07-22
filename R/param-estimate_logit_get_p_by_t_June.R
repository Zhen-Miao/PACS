## param_estimate_logit


## check if it returns error message
is.error <- function(x) inherits(x, "try-error")

## loss function -- original

#' Title Calculate Loss Function without Firth prior
#'
#' @param q_vec A vector of capturing rates
#' @param y_vec A vector of observed open (1) or close (0) state for each cell,
#'  must match the length of q
#' @param p_bg A vector of open probability
#'
#' @return Log loss value without Firth prior
#' @export
#'
loss_fun <- function(p_bg, q_vec, y_vec) {
  log_pq <- log(q_vec) + log(p_bg)
  log_1_pq <- log(1 - q_vec * p_bg)
  log_loss <- sum(y_vec * log_pq + (1 - y_vec) * log_1_pq)
  return(log_loss)
}




#' Title Gradient of loss function without Firth prior
#'
#' @param xdumm A design matrix, X with dummy variables
#' @param p_bg A value of open probability
#' @param q_vec A vector of capturing rates
#' @param y_vec A vector of observed open (1) or close (0) state for each cell,
#'  must match the length of q
#'
#' @return The first order derivative of the loss function without Firth prior
#' @export
loss_gradient <- function(xdumm, p_bg, q_vec, y_vec) {
  ## the first part
  y_m_p <- y_vec - p_bg
  grad_vec1 <- crossprod(xdumm, y_m_p)

  ## the second part
  fc <- (1 - q_vec) / (1 - p_bg * q_vec)
  y_p <- (1 - y_vec) * p_bg * fc
  grad_vec2 <- crossprod(xdumm, y_p)

  grad_vec <- as.vector(grad_vec1) + as.vector(grad_vec2)
  return(grad_vec)
}


#' Title Information matrix -- original without Firth prior
#'
#' @param xdumm A design matrix, X with dummy variables
#' @param p_bg A vector of open probability
#' @param q_vec A vector of capturing rates
#'
#' @return information matrix
#' @export
infor_mat <- function(xdumm, p_bg, q_vec) {
  wii <- (p_bg * q_vec * (1 - p_bg)^2) / (1 - p_bg * q_vec)
  wii_sqrt <- sqrt(as.vector(wii))
  wii_sqrt_X <- wii_sqrt * xdumm
  inf_m <- crossprod(wii_sqrt_X)

  return(inf_m)
}


#' Title Compute sqrt(Wii) * X
#'
#' @param xdumm Design matrix
#' @param p_bg open probability vector
#' @param q_vec Capturing probability vector
#'
#' @return sqrt(Wii) * X
#' @noRd
compute_wii_sqrt_X <- function(xdumm, p_bg, q_vec) {
  wii <- (p_bg * q_vec * (1 - p_bg)^2) / (1 - p_bg * q_vec)
  wii_sqrt <- sqrt(as.vector(wii))
  wii_sqrt_X <- wii_sqrt * xdumm

  return(wii_sqrt_X)
}


##
#' Title Gradient of penalized loss ( Invariant prior)
#'
#' @param xdumm A design matrix, X
#' @param p_bg A value of open probability
#' @param q_vec A vector of capturing probability
#' @param inf_mat Information matrix
#' @param wii_sqrt_X sqrt(wii) as a vector times the design matrix
#'
#' @return Gradient of penalized loss ( Invariant prior)
#' @export
loss_grad_pen <- function(xdumm, p_bg, q_vec, inf_mat, wii_sqrt_X) {
  kii <- (2 * p_bg^2 * q_vec - 3 * p_bg + 1) / (1 - p_bg * q_vec)
  k_vec <- as.vector(kii)

  ## get the inverse of the information matrix
  i_infor_m <- try(solve(inf_mat), silent = TRUE)
  if (is.error(i_infor_m)) {
    return(NA)
  }

  h_mat1 <- wii_sqrt_X %*% i_infor_m ## n * p
  h_mat_diag <- rowSums(h_mat1 * wii_sqrt_X)

  hk <- as.matrix(h_mat_diag * k_vec)

  return(1 / 2 * crossprod(xdumm, hk))
}



#' Title compute_information matrix tilda
#'
#' @param xdumm A design matrix, X
#' @param p_bg open probability vector
#' @param q_vec A vector of capturing probability
#' @param y_vec A vector of observed open (1) or close (0) state for each cell,
#'  must match the length of q
#'
#' @return information matrix tilda
#' @export
compute_infor_mat_tilda <- function(xdumm, p_bg, q_vec, y_vec) {
  wii <- ((-1 * p_bg^2 * q_vec^2) +
            q_vec * (2 * p_bg + y_vec - 1) - y_vec) * p_bg * (1 - p_bg) /
    ((1 - p_bg * q_vec)^2)
  wii_ <- -1 * as.vector(wii)
  wii_X <- wii_ * xdumm
  inf_m <- crossprod(xdumm, wii_X)
  return(inf_m)
}



#' Title Iteratively reweighted least squares function (main function)
#'
#' @param xdumm Design matrix
#' @param theta_estimated Estimated coefficient values
#' @param stop_criteria Resolution to stop, default = 1e-6
#' @param q_vec Capturing probability vector
#' @param y_vec Observed counts
#'
#' @return Converged coefficient estimate
#' @export
irls_iter <- function(
    y_vec, xdumm, theta_estimated,
    stop_criteria = 1e-6, q_vec
    ) {
  indi <- 1
  n_iter <- 1
  conv.stat <- 1
  ## 1 means converge, 2 means matrix singular, 3 means reaches max iter
  # y_vec = as.vector(y_vec)
  while (indi >= stop_criteria && n_iter <= 15) {
    ## compute some values
    x_times_theta <- xdumm %*% theta_estimated
    p_bg <- 1 - 1 / (exp(x_times_theta) + 1)

    wii_sqrt_X <- compute_wii_sqrt_X(xdumm, p_bg, q_vec = q_vec)
    inf_mat <- crossprod(wii_sqrt_X)

    if (determinant.matrix(inf_mat)$modulus[1] == 0) {
      conv.stat <- 2
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      return(theta_n_conv)
    }

    score_total <- loss_grad_pen(
      xdumm, p_bg, q_vec = q_vec, inf_mat, wii_sqrt_X
      ) +
      loss_gradient(xdumm, p_bg, q_vec, y_vec)

    if (!anyNA(score_total)) {
      inf_mat_inv <- try(solve(inf_mat), silent = TRUE)
    } else {
      conv.stat <- 2
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      return(theta_n_conv)
    }

    theta_update <- theta_estimated + inf_mat_inv %*% score_total
    indi <- sum((theta_update - theta_estimated)^2)
    theta_estimated <- theta_update
    n_iter <- n_iter + 1
  }

  if (indi >= stop_criteria) {
    conv.stat <- 3
  }

  theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
  return(theta_n_conv)
}


#' Title Iteratively reweighted least squares function (main function) for the
#'  null model
#'
#' @param y_vec Observed counts
#' @param xdumm Design matrix
#' @param theta_estimated Estimated coefficient values
#' @param hold_zero Which element to hold zero in the null model
#' @param stop_criteria Resolution to stop, default = 1e-6
#' @param q_vec Capturing probability vector
#'
#' @return Converged coefficient estimate for the null model
#' @export
irls_iter_null <- function(
    y_vec, xdumm, theta_estimated,
    hold_zero, stop_criteria = 1e-6, q_vec) {
  indi <- 1
  n_iter <- 1
  conv.stat <- 1
  ## 1 means converge, 2 means matrix singular, 3 means reaches max iter
  # y_vec = as.vector(y_vec)
  poi <- setdiff(c(seq_along(theta_estimated)), hold_zero)
  zero_augment <- rep(0, length(hold_zero))
  while (indi >= stop_criteria && n_iter <= 15) {
    ## compute some values
    x_times_theta <- xdumm %*% theta_estimated
    p_bg <- 1 - 1 / (exp(x_times_theta) + 1)

    wii_sqrt_X <- compute_wii_sqrt_X(xdumm, p_bg, q_vec = q_vec)
    inf_mat <- crossprod(wii_sqrt_X)

    if (determinant.matrix(inf_mat)$modulus[1] == 0) {
      conv.stat <- 2
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      return(theta_n_conv)
    }

    score_total <- loss_grad_pen(
      xdumm, p_bg, q_vec = q_vec, inf_mat, wii_sqrt_X
      ) +
      loss_gradient(xdumm, p_bg, q_vec, y_vec)

    if (!anyNA(score_total)) {
      inf_mat_inv <- try(solve(inf_mat[poi, poi]), silent = TRUE)
    } else {
      conv.stat <- 2
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      return(theta_n_conv)
    }

    updated_offset <- c(inf_mat_inv %*% score_total[poi], zero_augment)
    theta_update <- theta_estimated + updated_offset
    indi <- sum((theta_update - theta_estimated)^2)
    theta_estimated <- theta_update
    n_iter <- n_iter + 1
  }

  if (indi >= stop_criteria) {
    conv.stat <- 3
  }

  theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
  return(theta_n_conv)
}


#' Iteratively reweighted least squares function -- newton methods
#'
#' @param xdumm Design matrix
#' @param theta_estimated Estimated coefficient values
#' @param stop_criteria Resolution to stop, default = 1e-6
#' @param q_vec Capturing probability vector
#' @param y_vec Observed counts
#'
#' @return Converged coefficient estimate
#' @export
irls_iter_nt <- function(
    y_vec, xdumm, theta_estimated,
    stop_criteria = 1e-5, q_vec
    ) {
  indi <- 1
  n_iter <- 1
  conv.stat <- 1
  ## 1 means converge, 2 means matrix singular, 3 means reaches max iter
  # y_vec = as.vector(y_vec)
  loss_star <- vector(length = 15)
  while (indi >= stop_criteria && n_iter <= 15) {
    ## compute some values
    x_times_theta <- xdumm %*% theta_estimated
    p_bg <- 1 - 1 / (exp(x_times_theta) + 1)

    inf_mat_tilda <- compute_infor_mat_tilda(xdumm, p_bg, q_vec, y_vec)
    wii_sqrt_X <- compute_wii_sqrt_X(xdumm, p_bg, q_vec = q_vec)
    inf_mat <- crossprod(wii_sqrt_X)

    loss_star[n_iter] <- loss_fun_star(xdumm, p_bg, q_vec, y_vec)
    # print(loss_star[n_iter])
    ## note, we are still
    ## using the previous definition of w_ii_sqrt_X,
    ## this is because the other way of computing this value does not always
    ## result in positive values, which cannot
    ## proceed. If the two working weight both make sense,
    ## they should not be too different
    ## and thus the results should still converge with this combination

    inf_mat_tilda <- 0.5 * inf_mat_tilda + 0.5 * inf_mat
    # inf_mat_tilda = inf_mat

    if (determinant.matrix(inf_mat_tilda)$modulus[1] == 0) {
      conv.stat <- 2
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      return(theta_n_conv)
    }

    score_total <- loss_grad_pen(
      xdumm, p_bg, q_vec = q_vec, inf_mat_tilda, wii_sqrt_X
      ) +
      loss_gradient(xdumm, p_bg, q_vec, y_vec)

    if (!anyNA(score_total)) {
      inf_mat_inv <- try(solve(inf_mat_tilda), silent = TRUE)
    } else {
      conv.stat <- 2
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      return(theta_n_conv)
    }

    theta_update <- theta_estimated + inf_mat_inv %*% score_total
    # print(paste('theta update'))
    # print(theta_update)
    indi <- sum((theta_update - theta_estimated)^2)
    if (n_iter >= 2 &&
        !is.na(loss_star[n_iter]) &&
        loss_star[n_iter] != -Inf) {
      if (loss_star[n_iter] - loss_star[n_iter - 1] > 0) {
        theta_estimated <- theta_update
        if (loss_star[n_iter] - loss_star[n_iter - 1] < 1) {
          break
          ## stop updating if the change is not ver much,
          ## consider it reaches minimum
        }
      } else {
        # print('loss not increasing')
        conv.stat <- 4
        theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
        return(theta_n_conv)
        break
      }
    } else if (!is.na(loss_star[n_iter]) && loss_star[n_iter] != -Inf) {
      theta_estimated <- theta_update
    } else {
      conv.stat <- 4
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      break
    }
    n_iter <- n_iter + 1
  }

  if (n_iter > 15) {
    conv.stat <- 3
  }

  theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
  return(theta_n_conv)
}


#' Iteratively reweighted least squares function for null -- newton method
#'  for the null model
#'
#' @param y_vec Observed counts
#' @param xdumm Design matrix
#' @param theta_estimated Estimated coefficient values
#' @param hold_zero Which element to hold zero in the null model
#' @param stop_criteria Resolution to stop, default = 1e-6
#' @param q_vec Capturing probability vector
#'
#' @return Converged coefficient estimate for the null model
#' @export
irls_iter_nt_null <- function(
    y_vec, xdumm, theta_estimated, hold_zero,
    stop_criteria = 1e-5, q_vec
    ) {
  indi <- 1
  n_iter <- 1
  conv.stat <- 1 ## 1 means converge,
  ## 2 means matrix singular, 3 means reaches max iter
  # y_vec = as.vector(y_vec)
  poi <- setdiff(c(seq_along(theta_estimated)), hold_zero)
  zero_augment <- rep(0, length(hold_zero))
  loss_star <- vector(length = 15)
  while (indi >= stop_criteria && n_iter <= 15) {
    ## compute some values
    x_times_theta <- xdumm %*% theta_estimated
    p_bg <- 1 - 1 / (exp(x_times_theta) + 1)

    inf_mat_tilda <- compute_infor_mat_tilda(xdumm, p_bg, q_vec, y_vec)
    wii_sqrt_X <- compute_wii_sqrt_X(xdumm, p_bg, q_vec = q_vec)
    inf_mat <- crossprod(wii_sqrt_X)

    loss_star[n_iter] <- loss_fun_star(xdumm, p_bg, q_vec, y_vec)
    # print(loss_star[n_iter])
    ## note, we are still
    ## using the previous definition of w_ii_sqrt_X,
    ## this is because the other way of
    ## computing this value does not always result in positive values,
    ## which cannot
    ## proceed. If the two working weight both make sense,
    ## they should not be too different
    ## and thus the results should still coverge with this combination

    inf_mat_tilda <- 0.5 * inf_mat_tilda + 0.5 * inf_mat
    # inf_mat_tilda = inf_mat

    if (determinant.matrix(inf_mat_tilda)$modulus[1] == 0) {
      conv.stat <- 2
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      return(theta_n_conv)
    }

    score_total <- loss_grad_pen(
      xdumm, p_bg, q_vec = q_vec, inf_mat_tilda, wii_sqrt_X
      ) +
      loss_gradient(xdumm, p_bg, q_vec, y_vec)

    if (!anyNA(score_total)) {
      inf_mat_inv <- try(solve(inf_mat_tilda[poi, poi]), silent = TRUE)
    } else {
      conv.stat <- 2
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      return(theta_n_conv)
    }

    updated_offset <- c(inf_mat_inv %*% score_total[poi], zero_augment)
    theta_update <- theta_estimated + updated_offset

    # print(paste('theta update'))
    # print(theta_update)
    indi <- sum((theta_update - theta_estimated)^2)
    if (n_iter >= 2 &&
        !is.na(loss_star[n_iter]) &&
        loss_star[n_iter] != -Inf) {
      if (loss_star[n_iter] - loss_star[n_iter - 1] > 0) {
        theta_estimated <- theta_update
        if (loss_star[n_iter] - loss_star[n_iter - 1] < 1) {
          break
          ## stop updating if the change is not ver much,
          ## consider it reaches minimum
        }
      } else {
        # print('loss not increasing')
        conv.stat <- 4
        theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
        return(theta_n_conv)
        break
      }
    } else if (!is.na(loss_star[n_iter]) && loss_star[n_iter] != -Inf) {
      theta_estimated <- theta_update
    } else {
      conv.stat <- 4
      theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
      break
    }
    n_iter <- n_iter + 1
  }

  if (n_iter > 15) {
    conv.stat <- 3
  }

  theta_n_conv <- c(as.vector(theta_estimated), conv.stat)
  return(theta_n_conv)
}



#' Title Regularized log loss function
#'
#' @param xdumm Observed counts
#' @param p_bg Vector of open probability
#' @param q_vec Capturing probability vector
#' @param y_vec Observed counts
#'
#' @return Regularized log loss
#' @export
loss_fun_star <- function(xdumm, p_bg, q_vec, y_vec) {
  log_pq <- log(q_vec) + log(p_bg)
  log_1_pq <- log(1 - q_vec * p_bg)
  log_loss <- sum(y_vec * log_pq + (1 - y_vec) * log_1_pq)

  wii <- (p_bg * q_vec * (1 - p_bg)^2) / (1 - p_bg * q_vec)
  wii_sqrt <- sqrt(as.vector(wii))
  wii_sqrt_X <- wii_sqrt * xdumm
  inf_m <- crossprod(wii_sqrt_X)

  log_loss_star <- log_loss + 0.5 * log(det(inf_m))
  return(log_loss_star)
}



#' Estimate parameters for the full model
#'
#' @importFrom tictoc tic
#' @importFrom tictoc toc
#' @importFrom parallel mclapply
#' @import Matrix
#'
#' @param r_by_c Region by cell matrix
#' @param design_mat Desing matrix
#' @param par_initial Initialized parameters
#' @param cap_rate_vec A vector of capturing probability
#' @param mc.cores Number of multi-cores, default = 1
#'
#' @return Estimated coefficient for the full model
#' @export
estimate_parameters <- function(
    r_by_c,
    design_mat,
    par_initial,
    cap_rate_vec,
    mc.cores = 1) {
  ## for sparse matrix, firstly convert it to regular matrix
  if ("sparseMatrix" %in% is(r_by_c)) {
    r_by_c <- as.matrix(r_by_c)
  }
  n_features <- dim(r_by_c)[1]

  r_by_c <- t(r_by_c)

  r_by_c <- as.list(as.data.frame(r_by_c)) ## now it is cell by region matrix

  tictoc::tic("running parameter_estimation")

  theta_estimated_one_list <- mclapply(r_by_c,
    FUN = irls_iter,
    xdumm = design_mat, theta_estimated = par_initial,
    q_vec = cap_rate_vec, mc.cores = mc.cores
  )

  tictoc::toc()

  n_divide <- ceiling(n_features / 4)
  o1 <- do.call(cbind, theta_estimated_one_list[1:n_divide])
  o2 <- do.call(cbind,
                theta_estimated_one_list[(n_divide + 1):(n_divide * 2)])
  o3 <- do.call(cbind,
                theta_estimated_one_list[(n_divide * 2 + 1):(n_divide * 3)])
  o4 <- do.call(cbind,
                theta_estimated_one_list[(n_divide * 3 + 1):n_features])

  theta_summaries <- cbind(o1, o2, o3, o4)

  ## for those fail to converge, check newton method
  fail_converge <- theta_summaries[dim(theta_summaries)[1], ] != 1
  print(paste(sum(fail_converge), "regions fail to converge, try newton"))
  print(head(which(fail_converge)))

  fc_regions <- names(fail_converge[fail_converge])

  if (sum(fail_converge) == 0) {
    return(theta_summaries)
  } else {

    r_by_c <- r_by_c[fc_regions]

    par_initial <- -1 * par_initial
    tictoc::tic("newton method for estimation")

    theta_estimated_one_list_nt <- mclapply(r_by_c,
      FUN = irls_iter_nt,
      xdumm = design_mat, theta_estimated = par_initial,
      q_vec = cap_rate_vec, mc.cores = 1
    )
    tictoc::toc()

    if (length(theta_estimated_one_list_nt) == 1) {
      theta_summaries[, fc_regions] <- theta_estimated_one_list_nt[[1]]
    } else {
      theta_nt <- do.call(cbind, theta_estimated_one_list_nt)
      theta_summaries[, fc_regions] <- theta_nt[, fc_regions]
      fail_converge <- theta_nt[dim(theta_nt)[1], ] != 1
      print(paste(sum(fail_converge), "regions fail to converge aftere newton"))
    }
  }

  return(theta_summaries)
}



#' Estimate parameters for the null model
#' @importFrom tictoc tic
#' @importFrom tictoc toc
#' @importFrom parallel mclapply
#' @import Matrix
#'
#' @param r_by_c Region by cell matrix
#' @param design_mat Desing matrix
#' @param par_initial Initialized parameters
#' @param hold_zero Which parameter to hold zero during estimation
#' @param cap_rate_vec A vector of capturing probability
#' @param mc.cores Number of multi-cores, default = 1
#'
#' @return Estimated coefficient for the null model
#' @export
estimate_parameters_null <- function(
    r_by_c,
    design_mat,
    par_initial,
    hold_zero,
    cap_rate_vec,
    mc.cores = 1) {
  ## for sparse matrix, firstly convert it to regular matrix
  if ("sparseMatrix" %in% is(r_by_c)) {
    r_by_c <- as.matrix(r_by_c)
  }
  n_features <- dim(r_by_c)[1]

  r_by_c <- t(r_by_c)

  r_by_c <- as.list(as.data.frame(r_by_c)) ## now it is cell by region matrix

  tictoc::tic("parameter_estimation")
  ### to-do: add make clusters and so on to enable parallel computing


  theta_estimated_one_list <- mclapply(r_by_c,
    FUN = irls_iter_null,
    xdumm = design_mat, theta_estimated = par_initial,
    hold_zero = hold_zero,
    q_vec = cap_rate_vec, mc.cores = mc.cores
  )
  tictoc::toc()

  n_divide <- ceiling(n_features / 4)
  o1 <- do.call(cbind, theta_estimated_one_list[1:n_divide])
  o2 <- do.call(cbind,
                theta_estimated_one_list[(n_divide + 1):(n_divide * 2)])
  o3 <- do.call(cbind,
                theta_estimated_one_list[(n_divide * 2 + 1):(n_divide * 3)])
  o4 <- do.call(cbind,
                theta_estimated_one_list[(n_divide * 3 + 1):n_features])

  theta_summaries <- cbind(o1, o2, o3, o4)

  ## for those fail to converge, check newton method
  fail_converge <- theta_summaries[dim(theta_summaries)[1], ] != 1
  print(paste(sum(fail_converge), "regions fail to converge, try newton"))
  print(head(which(fail_converge)))

  fc_regions <- names(fail_converge[fail_converge])

  if (sum(fail_converge) == 0) {
    return(theta_summaries)
  } else {

    r_by_c <- r_by_c[fc_regions]

    par_initial <- -1 * par_initial
    tictoc::tic("newton method for estimation")

    theta_estimated_one_list_nt <- mclapply(r_by_c,
      FUN = irls_iter_nt_null,
      xdumm = design_mat, theta_estimated = par_initial,
      hold_zero = hold_zero,
      q_vec = cap_rate_vec, mc.cores = 1
    )
    tictoc::toc()

    if (length(theta_estimated_one_list_nt) == 1) {
      theta_summaries[, fc_regions] <- theta_estimated_one_list_nt[[1]]
    } else {
      theta_nt <- do.call(cbind, theta_estimated_one_list_nt)
      theta_summaries[, fc_regions] <- theta_nt[, fc_regions]
      fail_converge <- theta_nt[dim(theta_nt)[1], ] != 1
      print(paste(sum(fail_converge),
                  "regions fail to converge aftere newton"))
    }
  }

  return(theta_summaries)
}
