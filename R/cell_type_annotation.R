## Functions for cell type label annotation

#' Estimate cell type label for new data
#'
#' @importFrom Rfast colsums
#' @import Matrix
#'
#' @param r_by_t Region-by cell type matrix generated from a annotated dataset
#' @param in_r_by_c Input (unannotated) region by cell matrix
#' @param alpha Weight for negative peaks, default = 1
#'
#' @return A matrix of cell types by cell matrix, with elements representing
#'  probability of being in that cell type
#' @export
estimate_label_default <- function(
    r_by_t,
    in_r_by_c,
    alpha = 1) {
  ## build prediction matrix

  ## row cell, column cell types
  esti_m5 <- matrix(nrow = dim(in_r_by_c)[2], ncol = dim(r_by_t)[2])
  rownames(esti_m5) <- colnames(in_r_by_c)
  colnames(esti_m5) <- colnames(r_by_t)

  if (inherits(in_r_by_c, "Matrix")) {
    if (any(in_r_by_c@x != 0 & in_r_by_c@x != 1)) {
      stop("the matrix should be a binary matrix!")
    }
  } else {
    if (any(in_r_by_c != 0 & in_r_by_c != 1)) {
      stop("the matrix should be a binary matrix!")
    }
  }
  n_reads <- Matrix::colSums(in_r_by_c)
  sum_prob <- colSums(r_by_t)
  capturing_rate_mat <- t(outer(n_reads, sum_prob, "/"))
  capturing_rate_mat[capturing_rate_mat > 0.9995] <- 0.9995
  capturing_rate_mat[capturing_rate_mat < 0.00005] <- 0.00005

  if (as.numeric(dim(in_r_by_c)[1]) * as.numeric(dim(in_r_by_c)[2])
      > (2^31 - 1)) {
    ## need to split
    split_by <- ceiling((as.numeric(dim(in_r_by_c)[1]) *
                           as.numeric(dim(in_r_by_c)[2])) / (2^31 - 1))
    dim2 <- ceiling(dim(in_r_by_c)[2] / split_by)
    for (i in c(1:split_by)) {
      in_r_by_c_mat <- as.matrix(
        in_r_by_c[, ((i - 1) * dim2 + 1):(min(i * dim2, dim(in_r_by_c)[2]))]
        )

      ## predict cell types
      for (i_cell in seq_len(dim(in_r_by_c_mat)[2])) {
        pqbyt <- r_by_t * capturing_rate_mat[, i_cell + (i - 1) * dim2]
        lg_pqbyt <- log(pqbyt)
        lg_pqbyt_q <- log1p(-1 * pqbyt)
        # unparallelized
        xvec <- in_r_by_c_mat[, i_cell]
        esti_m5[i_cell + (i - 1) * dim2, ] <-
          Rfast::colsums(lg_pqbyt * xvec + alpha * lg_pqbyt_q * (1 - xvec))
      }
      saveRDS(esti_m5, "esti_m5_intermediate_ave_high q_12.rds")
    }
  } else {
    in_r_by_c_mat <- as.matrix(in_r_by_c)

    ## predict cell types

    # unparallelized
    for (i_cell in seq_len(dim(esti_m5)[1])) {
      pqbyt <- r_by_t * capturing_rate_mat[, i_cell]
      lg_pqbyt <- log(pqbyt)
      lg_pqbyt_q <- log1p(-1 * pqbyt)
      # unparallelized
      xvec <- in_r_by_c_mat[, i_cell]
      esti_m5[i_cell, ] <- Rfast::colsums(lg_pqbyt * xvec +
                                            lg_pqbyt_q * (1 - xvec))

      # esti_m5[i_cell,] =
      # apply(lg_pqbyt, 2, function(x) sum(x[in_r_by_c_mat[,i_cell] > 0])) +
      #   alpha*  apply(lg_pqbyt_q, 2, function(x)
      # sum(x[in_r_by_c_mat[,i_cell] == 0]))
    }
    return(esti_m5)
  }
}

#' Estimate cell type label for new data using all peaks as input
#'
#' @importFrom Rfast colsums
#' @import Matrix
#' @param r_by_t Region-by cell type matrix generated from a annotated dataset
#' @param in_r_by_c Input (unannotated) region by cell matrix
#' @param alpha Weight for negative peaks, default = 1
#' @param pks_sel selected peaks as informative for cell type label prediction
#'
#' @return A matrix of cell types by cell matrix, with elements representing
#'  probability of being in that cell type
#' @export
estimate_label_no_cap_rate_all_pks <- function(
    r_by_t,
    in_r_by_c,
    pks_sel,
    alpha = 1) {
  ## build prediction matrix

  ## row cell, column cell types
  esti_m5 <- matrix(nrow = dim(in_r_by_c)[2], ncol = dim(r_by_t)[2])
  rownames(esti_m5) <- colnames(in_r_by_c)
  colnames(esti_m5) <- colnames(r_by_t)

  if (inherits(in_r_by_c, "Matrix")) {
    if (any(in_r_by_c@x != 0 & in_r_by_c@x != 1)) {
      stop("the matrix should be a binary matrix!")
    }
  } else {
    if (any(in_r_by_c != 0 & in_r_by_c != 1)) {
      stop("the matrix should be a binary matrix!")
    }
  }
  n_reads <- Matrix::colSums(in_r_by_c)
  sum_prob <- colSums(r_by_t)
  capturing_rate_mat <- t(outer(n_reads, sum_prob, "/"))
  capturing_rate_mat[capturing_rate_mat > 0.9995] <- 0.9995
  capturing_rate_mat[capturing_rate_mat < 0.00005] <- 0.00005

  ## filter peaks to keep only relevant ones
  r_by_t <- r_by_t[pks_sel, ]
  in_r_by_c <- in_r_by_c[pks_sel, ]

  if (as.numeric(dim(in_r_by_c)[1]) * as.numeric(dim(in_r_by_c)[2]) > (2^31 - 1)) {
    ## need to split
    split_by <- ceiling((as.numeric(dim(in_r_by_c)[1]) * as.numeric(dim(in_r_by_c)[2])) / (2^31 - 1))
    dim2 <- ceiling(dim(in_r_by_c)[2] / split_by)
    for (i in c(1:split_by)) {
      in_r_by_c_mat <- as.matrix(in_r_by_c[, ((i - 1) * dim2 + 1):(min(i * dim2, dim(in_r_by_c)[2]))])
      # in_capturing_rate = capturing_rate[((i-1)*dim2+1):(min(i*dim2,dim(in_r_by_c)[2]))]

      ## predict cell types
      for (i_cell in seq_len(dim(in_r_by_c_mat)[2])) {
        pqbyt <- r_by_t * capturing_rate_mat[, i_cell + (i - 1) * dim2]
        lg_pqbyt <- log(pqbyt)
        lg_pqbyt_q <- log1p(-1 * pqbyt)
        # unparallelized
        xvec <- in_r_by_c_mat[, i_cell]
        esti_m5[i_cell + (i - 1) * dim2, ] <- Rfast::colsums(lg_pqbyt * xvec + alpha * lg_pqbyt_q * (1 - xvec))
      }
      saveRDS(esti_m5, "esti_m5_intermediate_ave_high q_12.rds")
    }
  } else {
    in_r_by_c_mat <- as.matrix(in_r_by_c)

    ## predict cell types

    # unparallelized
    for (i_cell in seq_len(dim(esti_m5)[1])) {
      pqbyt <- r_by_t * capturing_rate_mat[, i_cell]
      lg_pqbyt <- log(pqbyt)
      lg_pqbyt_q <- log1p(-1 * pqbyt)
      # unparallelized
      xvec <- in_r_by_c_mat[, i_cell]
      esti_m5[i_cell, ] <- Rfast::colsums(lg_pqbyt * xvec + lg_pqbyt_q * (1 - xvec))

      # esti_m5[i_cell,] = apply(lg_pqbyt, 2, function(x) sum(x[in_r_by_c_mat[,i_cell] > 0])) +
      #   alpha*  apply(lg_pqbyt_q, 2, function(x) sum(x[in_r_by_c_mat[,i_cell] == 0]))
    }
    return(esti_m5)
  }
}
