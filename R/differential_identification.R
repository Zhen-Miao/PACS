### differential_identification

## remember to load the LRT functions
if (!exists("CCT_internal")) {
  stop("please load the Cauchy combination test function first!")
}

#' Title Likely hood ratio test between a cell type and background cell type
#'
#' @param ir iteration indicator
#' @param q_vec A vector of capturing rate
#' @param p_by_c Peak by cell matrix
#' @param x_one Design matrix for the full model
#' @param par_initial_one Initial parameter point for full model
#' @param x_null Design matrix for the null model
#' @param par_initial_null Initial parameter point for null model
#' @param stop_criteria Criteria for determining convergence
#' @param n_feature Number of features (peaks or bins)
#' @param n_group Number of groups
#' @param group_set The set of groups
#'
#' @return Differential P values for each feature and its convergence indicator
#' @export
#'
#' @examples
LRT_one_vs_bg <- function(ir,
                          q_vec, p_by_c,
                          x_one, par_initial_one,
                          x_null, par_initial_null,
                          stop_criteria = 1e-6,
                          n_feature,
                          n_group,
                          group_set) {

  ## initiation
  conv.stat <- 1 ## 1 means converge, 2 means matrix singular, 3 means reaches max iter
  ## 4 means singular during null , 5 means max iter during null
  LRT_p_mat <- vector(length = 2)
  names(LRT_p_mat) <- c("p_value", "conv_stat")


  ## get y vector
  ir_c <- as.character(ir)
  y_vec <- p_by_c[ir, ]

  theta_estimated_one_list <- rwls_iter(
    x_one, par_initial_one,
    stop_criteria, q_vec, y_vec
  )
  theta_estimated_one <- theta_estimated_one_list$thata_final
  theta_one_conv <- theta_estimated_one_list$conv.stat

  if (theta_one_conv != 1) {
    theta_estimated_one_list <- rwls_iter_nt(x_one, par_initial_one, stop_criteria, q_vec, y_vec)
    theta_estimated_one <- theta_estimated_one_list$thata_final
    theta_one_conv <- theta_estimated_one_list$conv.stat
    if (theta_one_conv != 1) {
      LRT_p_mat["conv_stat"] <- theta_one_conv
      LRT_p_mat["p_value"] <- NA
      return(LRT_p_mat)
    }
  }

  theta_estimated_null_list <- rwls_iter(
    x_null, par_initial_null,
    stop_criteria, q_vec, y_vec
  )
  theta_estimated_null <- theta_estimated_null_list$thata_final
  theta_null_conv <- theta_estimated_null_list$conv.stat

  if (theta_null_conv != 1) {
    theta_estimated_null_list <- rwls_iter_nt(x_null, par_initial_null, stop_criteria, q_vec, y_vec)
    theta_estimated_null <- theta_estimated_null_list$thata_final
    theta_null_conv <- theta_estimated_null_list$conv.stat

    if (theta_null_conv != 1) {
      LRT_p_mat["conv_stat"] <- theta_one_conv
      LRT_p_mat["p_value"] <- NA
      return(LRT_p_mat)
    }
  }

  theta_estimated_null_3 <- append(theta_estimated_null, 0)
  LR_value <- 2 * (loss_fun_star(x_one, theta_estimated_one, q_vec, y_vec) -
    loss_fun_star(x_one, theta_estimated_null_3, q_vec, y_vec))

  LRT_p_value <- pchisq(LR_value, df = 1, lower.tail = FALSE)
  LRT_p_mat["conv_stat"] <- 1
  LRT_p_mat["p_value"] <- LRT_p_value

  return(LRT_p_mat)
}


LRT_one_vs_ot <- function(ir,
                          target,
                          background_set,
                          q_vec_f,
                          p_by_c,
                          x_one_f,
                          par_initial_one,
                          x_null_f,
                          par_initial_null,
                          stop_criteria = 1e-6,
                          n_feature,
                          n_group,
                          group_set,
                          thin_by_group,
                          combination_test_method = "cct",
                          verbose = F) {

  ## check input
  if (!(combination_test_method %in% c("cct", "max", "half", "min"))) {
    stop("please specify combine test method as one of c('cct', 'max','half','min')")
  }

  ## get the y vector
  y_vec_f <- p_by_c[ir, ]

  ## get some quantities
  n_bg <- length(background_set)

  ## initiation
  conv.stat <- 1 ## 1 means converge, 2 means matrix singular,
  ## 3 means reaches max iter
  ## 4 means singular during null , 5 means max iter during null
  LRT_p_cct <- vector(length = 2)
  names(LRT_p_cct) <- c("p_value", "conv_stat")

  ## intermediate check
  LRT_p_interm <- matrix(nrow = 2, ncol = n_bg)
  colnames(LRT_p_interm) <- background_set
  rownames(LRT_p_interm) <- c("p_value", "conv_stat")

  LRT_p_interm["conv_stat", ] <- 1 ## default assume all converge


  for (jj in 1:n_bg) {
    y_vec <- y_vec_f[thin_by_group[, jj]]
    x_one <- x_one_f[thin_by_group[, jj], ]
    q_vec <- q_vec_f[thin_by_group[, jj]]
    x_null <- x_null_f[thin_by_group[, jj], ]

    theta_estimated_one_list <- rwls_iter(
      x_one, par_initial_one,
      stop_criteria, q_vec, y_vec
    )
    theta_estimated_one <- theta_estimated_one_list$thata_final
    theta_one_conv <- theta_estimated_one_list$conv.stat

    if (theta_one_conv != 1) {
      theta_estimated_one_list <- rwls_iter_nt(x_one, par_initial_one, stop_criteria, q_vec, y_vec)
      theta_estimated_one <- theta_estimated_one_list$thata_final
      theta_one_conv <- theta_estimated_one_list$conv.stat
      if (theta_one_conv != 1) {
        LRT_p_interm["conv_stat", jj] <- theta_one_conv
        LRT_p_interm["p_value", jj] <- NA
        next
      }
    }

    theta_estimated_null_list <- rwls_iter(
      x_null, par_initial_null,
      stop_criteria, q_vec, y_vec
    )
    theta_estimated_null <- theta_estimated_null_list$thata_final
    theta_null_conv <- theta_estimated_null_list$conv.stat

    if (theta_null_conv != 1) {
      theta_estimated_null_list <- rwls_iter_nt(x_null, par_initial_null, stop_criteria, q_vec, y_vec)
      theta_estimated_null <- theta_estimated_null_list$thata_final
      theta_null_conv <- theta_estimated_null_list$conv.stat

      if (theta_null_conv != 1) {
        LRT_p_interm["conv_stat", jj] <- theta_one_conv
        LRT_p_interm["p_value", jj] <- NA
        next
      }
    }

    theta_estimated_null_3 <- append(theta_estimated_null, 0)
    LR_value <- 2 * (loss_fun_star(x_one, theta_estimated_one, q_vec, y_vec) -
      loss_fun_star(x_one, theta_estimated_null_3, q_vec, y_vec))

    LRT_p_value <- pchisq(LR_value, df = 1, lower.tail = FALSE)

    # LRT_p_matrix['conv_stat',ir] <- 1
    LRT_p_interm["p_value", jj] <- LRT_p_value
  }


  if (anyNA(LRT_p_interm["p_value", ])) {
    LRT_p_cct["p_value"] <- NA
    LRT_p_cct["conv_stat"] <- 10
  } else if (combination_test_method == "cct") {
    LRT_p_cct["p_value"] <- CCT_internal(LRT_p_interm["p_value", ])
    LRT_p_cct["conv_stat"] <- ifelse(test = all(LRT_p_interm["conv_stat", ]) == 1,
      yes = 1, no = 20
    )
  } else if (combination_test_method == "max") {
    LRT_p_cct["p_value"] <- max(LRT_p_interm["p_value", ])
    LRT_p_cct["conv_stat"] <- ifelse(test = all(LRT_p_interm["conv_stat", ]) == 1,
      yes = 1, no = 20
    )
  } else if (combination_test_method == "min") {
    LRT_p_cct["p_value"] <- min(LRT_p_interm["p_value", ])
    LRT_p_cct["conv_stat"] <- ifelse(test = all(LRT_p_interm["conv_stat", ]) == 1,
      yes = 1, no = 20
    )
  } else if (combination_test_method == "half") {
    LRT_p_cct["p_value"] <- quantile(LRT_p_interm["p_value", ], p = 0.5, type = 1)
    LRT_p_cct["conv_stat"] <- ifelse(test = all(LRT_p_interm["conv_stat", ]) == 1,
      yes = 1, no = 20
    )
  }

  return(LRT_p_cct)
  # if(ir == 20000 |  ir == 40000){
  #   saveRDS(object = LRT_p_matrix,
  #           file = paste('DAR_LR kidney PT alternative approach', ir, '.rds', sep = '_'))
  #
  # }
}



loss_firth_from_wii <- function(wii_sqrt, xdumm){
  wii_sqrt_X <- wii_sqrt * xdumm
  inf_mat <- crossprod(wii_sqrt_X)
  log_loss_star = 0.5 * determinant.matrix(inf_mat)$modulus[1]
  return(log_loss_star)
}



compare_models <- function(x_full, theta_estimated_full,
                           x_null, theta_estimated_null,
                           q_vec, c_by_r,df_test,mc.cores = 1
){

  if("sparseMatrix" %in% is(c_by_r)){
    c_by_r = as.matrix(c_by_r)
  }

  ## get number of regions (peaks)
  n_regions = dim(theta_estimated_full)[2]
  if(n_regions != dim(c_by_r)[2]){
    stop("please make sure the input matrix and the parameters match dimensions")
  }

  ## df
  if(is.vector(theta_estimated_null)){
    theta_estimated_null = matrix(theta_estimated_null, nrow = 1)
  }

  ###############
  ### for full matrix
  ###############

  ## vectorize -- calculate for all regions
  x_times_theta = x_full %*% theta_estimated_full ## n by r (number of features)
  p_bg <- 1 - 1 / (exp(x_times_theta) + 1) ## n by r (number of features)

  ## calculate loss without Firth prior
  log_pq <- log(q_vec) + log(p_bg) ## n by r (number of features)
  log_1_pq <- log(1 - p_bg * q_vec) ## n by r (number of features)
  log_loss = colSums(c_by_r * log_pq + (1 - c_by_r) * log_1_pq) ## vector, r

  ## calculate W for GLM
  wii <- (p_bg * (1 - p_bg) * (1 - p_bg) * q_vec) / (1 - p_bg * q_vec) ## n by r
  wii_sqrt <- sqrt(wii)  ## n by r

  ###############
  ### for null matrix
  ###############

  ## vectorize -- calculate for all regions
  x_times_theta_null = x_null %*% theta_estimated_null ## n by r (number of features)
  p_bg_null <- 1 - 1 / (exp(x_times_theta_null) + 1) ## n by r (number of features)

  ## calculate loss without Firth prior
  log_pq_null <- log(q_vec) + log(p_bg_null) ## n by r (number of features)
  log_1_pq_null <- log(1 - p_bg_null * q_vec) ## n by r (number of features)
  log_loss_null = colSums(c_by_r * log_pq_null + (1 - c_by_r) * log_1_pq_null) ## vector, r

  ## calculate W for GLM
  wii_null <- (p_bg_null * (1 - p_bg_null) * (1 - p_bg_null) * q_vec) / (1 - p_bg_null * q_vec) ## n by r
  wii_sqrt_null <- sqrt(wii_null)  ## n by r

  ## for each column, do this
  tic('computing firth loss start')
  log_loss_firth_full <- mclapply(as.data.frame(wii_sqrt),loss_firth_from_wii,
                                  xdumm = x_full, mc.cores = mc.cores)
  log_loss_firth_null <- mclapply(as.data.frame(wii_sqrt_null),loss_firth_from_wii,
                                  xdumm = x_null, mc.cores = mc.cores)
  toc()

  log_loss_star_full <- unlist(log_loss_firth_full) + log_loss ## vector length r
  log_loss_star_null <- unlist(log_loss_firth_null) + log_loss_null ## vector length r

  test_stat = 2 * (log_loss_star_full - log_loss_star_null)
  LRT_p_value <- pchisq(test_stat, df = df_test, lower.tail = FALSE)

  return(LRT_p_value)

}



