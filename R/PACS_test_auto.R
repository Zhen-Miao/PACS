## wrapper for covariate matrix


#' Profile likelihood ratio test with PACS
#'
#' @importFrom tictoc tic
#' @importFrom tictoc toc
#' @importFrom parallel mclapply
#'
#' @param covariate_meta.data A data.frame with columns representing the
#'   covariates and rows representing cells
#' @param formula_full A formula object representing the full model. For
#'   example, ~ cell_type + batch
#' @param formula_null A formula object representing the null model. For
#'   example, ~ batch
#' @param pic_matrix The input region-by-cell PIC matrix
#' @param n_peaks_per_round The number of peaks tested in each round.
#'   Default = NULL, which will be determined by the maximum allowable size of
#'   matrix in R.
#' @param T_proportion_cutoff To determine the maximum value of T, we set
#'   a criteria based on the proportion of reads that are >= T. This parameter
#'   set this proportion. Default = 0.2
#' @param cap_rates A vector of capturing probability for each cell
#' @param par_initial_null Initialized values of estimated parameters for the
#'   null model, we do not
#'   need to specify unless there are reasons to do so. Default = NULL
#' @param par_initial_full Initialized values of estimated parameters for the
#'   null model, we do not
#'   need to specify unless there are reasons to do so. Default = NULL
#' @param n_cores number of cores for multi-core computation
#'
#' @return A list of two elements, pacs_converged is a vector of length
#'   2*n_peaks representing the convergence status of the peak
#'   in the null and full model
#'   and pacs_p_val is a vector of length n_peaks representing the p values for
#'   each peak.
#' @export
#'
pacs_test_sparse <- function(covariate_meta.data, formula_full,
                             formula_null, pic_matrix,
                             n_peaks_per_round = NULL,
                             T_proportion_cutoff = 0.2,
                             cap_rates, par_initial_null = NULL,
                             par_initial_full = NULL, n_cores = 1,
                             verbose = TRUE) {
  ## check the data
  n_cell <- ncol(pic_matrix)
  n_peaks <- nrow(pic_matrix)
  p_names <- rownames(pic_matrix)

  if (nrow(covariate_meta.data) != n_cell) {
    stop("number of cells do not match between meta.data and data matrix")
  }

  if (length(cap_rates) != n_cell) {
    stop("number of cells do not match between cap_rates and data matrix")
  }

  if (is.null(rownames(pic_matrix))) {
    print("peak names not supplied, set to f_1 to f_n")
    rownames(pic_matrix) <- paste("f", 1:n_peaks, sep = "_")
  }

  if (is.null(n_peaks_per_round)) {
    n_peaks_per_round <- min(floor(2^30 / n_cell), n_peaks)
  }

  ########
  require("Matrix", quietly = TRUE)

  if (inherits(pic_matrix, "Matrix")) {
    ## counts that are >= 2
    pic_matrix_2 <- pic_matrix
    pic_matrix_2@x[pic_matrix_2@x == 1] <- 0

    pic_matrix_2 <- drop0(pic_matrix_2)
    pic_matrix_2@x <- rep(1, length = length(pic_matrix_2@x))

    ## binarized matrix
    pic_matrixbin <- pic_matrix
    pic_matrixbin@x <- rep(1, length = length(pic_matrixbin@x))
  } else {
    ## counts that are >= 2
    pic_matrix_2 <- Matrix(pic_matrix, sparse = TRUE)
    pic_matrix_2@x[pic_matrix_2@x == 1] <- 0

    pic_matrix_2 <- drop0(pic_matrix_2)
    pic_matrix_2@x <- rep(1, length = length(pic_matrix_2@x))

    ## binarized matrix
    pic_matrixbin <- Matrix(pic_matrix, sparse = TRUE)
    pic_matrixbin@x <- rep(1, length = length(pic_matrixbin@x))
  }

  ## calculate the proportion of counts >= 2
  rs <- rowSums(pic_matrixbin)
  rs2 <- rowSums(pic_matrix_2)

  p_2 <- rs2 / rs ## proportion of counts >= 2
  p_2[is.na(p_2)] <- 0
  # quantile(p_2, na.rm = T)

  n_p_2 <- sum(p_2 >= T_proportion_cutoff)
  n_p_b <- sum(p_2 < T_proportion_cutoff)

  if (verbose) {
    print(paste(n_p_2, "peaks consider cumulative logit models", sep = " "))
    print(paste(n_p_b, "peaks consider logit models", sep = " "))
  }

  rm(pic_matrix_2)
  gc(verbose = F)

  ## select features based on the proportion of 2
  f_sel <- names(p_2)[p_2 >= T_proportion_cutoff]
  f_b_sel <- names(p_2)[p_2 < T_proportion_cutoff]

  if (n_p_2 >= 1) {
    pic_matrix_cumu <- pic_matrix[f_sel, , drop = FALSE]

    n_iters <- ceiling(n_p_2 / n_peaks_per_round)
    p_cumu <- rep(list(), length = n_iters)

    for (jj in 1:n_iters) {
      peak_start <- (jj - 1) * n_peaks_per_round + 1
      peak_end <- min(n_p_2, jj * n_peaks_per_round)

      pic_dense <- as.matrix(pic_matrix_cumu[peak_start:peak_end, ])

      tic("pacs computing cumulative logit part")

      p_cumu[[jj]] <- pacs_test_cumu(
        covariate_meta.data = covariate_meta.data, max_T = 2,
        formula_full = formula_full,
        formula_null = formula_null,
        pic_matrix = pic_dense,
        cap_rates = cap_rates, n_cores = n_cores
      )
      toc()
      rm(pic_dense)
      gc(verbose = FALSE)
    }


    ## we do not need the quantitative part of the matrix any more
    rm(pic_matrix_cumu)
    rm(pic_matrix)
    gc(verbose = FALSE)
  }


  ## compute the binary part
  if (n_p_b >= 1) {
    pic_matrix_b <- pic_matrixbin[f_b_sel, , drop = F]

    n_iters_b <- ceiling(n_p_b / n_peaks_per_round)
    p_logit <- rep(list(), length = n_iters_b)

    for (jj in 1:n_iters_b) {
      peak_start <- (jj - 1) * n_peaks_per_round + 1
      peak_end <- min(n_p_b, jj * n_peaks_per_round)

      pic_dense <- as.matrix(pic_matrix_b[peak_start:peak_end, ])

      tic("pacs computing cumulative logit part")

      p_logit[[jj]] <- pacs_test_logit(
        covariate_meta.data = covariate_meta.data,
        formula_full = formula_full,
        formula_null = formula_null,
        pic_matrix = pic_dense,
        cap_rates = cap_rates, n_cores = n_cores
      )
      toc()
      rm(pic_dense)
      gc(verbose = F)
    }
  }


  ## organize the p values as well as convergence status
  if (n_p_2 >= 1 & n_p_b >= 1) {
    ## cumulative part
    if (n_iters > 1) {
      ## p values
      p_cumu_list <- sapply(p_cumu, function(x) x$pacs_p_val)
      p_val_cumu <- unlist(p_cumu_list)

      ## convergence status
      conv_cumu_list <- sapply(p_cumu, function(x)
        matrix(x$pacs_converged, ncol = 2))
      conv_cumu <- do.call(rbind, conv_cumu_list)
    } else {
      ## p values
      p_cumu_list <- sapply(p_cumu, function(x) x$pacs_p_val)
      p_val_cumu <- p_cumu_list[, 1, drop = T]

      ## convergence status
      conv_cumu_list <- sapply(p_cumu, function(x) x$pacs_converged)
      conv_cumu <- matrix(conv_cumu_list, ncol = 2)
      rownames(conv_cumu) <- names(p_val_cumu)
    }

    ## logit part
    if (n_iters_b > 1) {
      ## p values
      p_logit_list <- sapply(p_logit, function(x) x$pacs_p_val)
      p_val_logit <- unlist(p_logit_list)

      ## convergence status
      conv_logit_list <-
        sapply(p_logit, function(x) matrix(x$pacs_converged, ncol = 2))
      conver_logit <- do.call(rbind, conv_logit_list)
    } else {
      ## p values
      p_logit_list <- sapply(p_logit, function(x) x$pacs_p_val)
      p_val_logit <- p_logit_list[, 1, drop = TRUE]

      ## convergence status
      conv_logit_list <- sapply(p_logit, function(x) x$pacs_converged)
      conver_logit <- matrix(conv_logit_list, ncol = 2)
      rownames(conver_logit) <- names(p_val_logit)
    }

    p_val <- c(p_val_cumu, p_val_logit)[p_names]
    convergence <- rbind(conv_cumu, conver_logit)[p_names, ]
  } else if (n_p_2 == 0) {
    if (n_iters_b > 1) {
      ## p values
      p_logit_list <- sapply(p_logit, function(x) x$pacs_p_val)
      p_val_logit <- unlist(p_logit_list)
      p_val <- p_val_logit[p_names]

      ## convergence status
      conv_logit_list <-
        sapply(p_logit, function(x) matrix(x$pacs_converged, ncol = 2))
      conver_logit <- do.call(rbind, conv_logit_list)
      convergence <- conver_logit[p_names, ]
    } else {
      ## p values
      p_logit_list <- sapply(p_logit, function(x) x$pacs_p_val)
      p_val_logit <- p_logit_list[, 1, drop = TRUE]
      p_val <- p_val_logit[p_names]

      ## convergence status
      conv_logit_list <-
        sapply(p_logit, function(x) matrix(x$pacs_converged, ncol = 2))
      conver_logit <- matrix(conv_logit_list, ncol = 2)
      rownames(conver_logit) <- names(p_val_logit)
      convergence <- conver_logit[p_names, ]
    }
  } else if (n_p_b == 0) {
    if (n_iters > 1) {
      ## p values
      p_cumu_list <- sapply(p_cumu, function(x) x$pacs_p_val)
      p_val_cumu <- unlist(p_cumu_list)
      p_val <- p_val_cumu[p_names]

      ## convergence status
      conv_cumu_list <-
        sapply(p_cumu, function(x) matrix(x$pacs_converged, ncol = 2))
      conv_cumu <- do.call(rbind, conv_cumu_list)
      convergence <- conv_cumu[p_names, ]
    } else {
      ## p values
      p_cumu_list <- sapply(p_cumu, function(x) x$pacs_p_val)
      p_val_cumu <- p_cumu_list[, 1, drop = TRUE]
      p_val <- p_val_cumu[p_names]

      ## convergence status
      conv_cumu_list <- sapply(p_cumu, function(x) x$pacs_converged)
      conv_cumu <- matrix(conv_cumu_list, ncol = 2)
      rownames(conv_cumu) <- names(p_val_cumu)
      convergence <- conv_cumu[p_names, ]
    }
  }


  return(list(pacs_converged = convergence, pacs_p_val = p_val))
}
