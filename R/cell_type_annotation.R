## Functions for cell type label annotation

#' Estimate cell type label for new data
#'
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
  ## input check
  if (inherits(in_r_by_c, "Matrix")) {
    if (any(in_r_by_c@x != 0 & in_r_by_c@x != 1)) {
      stop("the matrix should be a binary matrix!")
    }
  } else {
    if (any(in_r_by_c != 0 & in_r_by_c != 1)) {
      stop("the matrix should be a binary matrix!")
    }
  }

  # Ensure r_by_t is a matrix
  r_by_t <- as.matrix(r_by_t)

  # Define the dimensions for easier reference
  n_cells <- ncol(in_r_by_c)
  cell_type_names <- colnames(r_by_t)

  # Pre-calculate sums to save computational time
  n_reads <- colSums(in_r_by_c)
  sum_prob <- colSums(r_by_t)
  capturing_rate_mat <- t(outer(n_reads, sum_prob, "/"))
  capturing_rate_mat[capturing_rate_mat > 0.9995] <- 0.9995
  capturing_rate_mat[capturing_rate_mat < 0.00005] <- 0.00005

  ## nested function to estimate the cell types
  cell_by_type_prob <- function(i_cell) {
    pqbyt <- r_by_t * capturing_rate_mat[, i_cell]
    lg_pqbyt <- log(pqbyt)
    lg_pqbyt_q <- log1p(-1 * pqbyt)
    xvec <- in_r_by_c_mat[, i_cell]
    esti_ct <- colSums(lg_pqbyt * xvec + lg_pqbyt_q * (1 - xvec) * alpha)
    return(esti_ct)
  }

  if (as.numeric(dim(in_r_by_c)[1]) * as.numeric(n_cells) > (2^31 - 1)) {
    ## need to split
    split_by <- ceiling((as.numeric(dim(in_r_by_c)[1]) *
      as.numeric(n_cells)) / (2^31 - 1))
    dim2 <- ceiling(dim(in_r_by_c)[2] / split_by)

    esti_m <- rep(list(), length = split_by)

    for (i in c(1:split_by)) {
      in_r_by_c_mat <- as.matrix(
        in_r_by_c[, ((i - 1) * dim2 + 1):(min(i * dim2, n_cells))]
      )

      ## predict cell types
      esti_list <- lapply(seq_len(dim(in_r_by_c_mat)[2]), cell_by_type_prob)
      esti_m[[i]] <- do.call(rbind(esti_list))
      rownames(esti_m[[i]]) <- colnames(in_r_by_c)
      colnames(esti_m[[i]]) <- cell_type_names
    }

    ## aggregate the results
    esti_mat <- do.call(rbind, esti_m)
  } else {
    in_r_by_c_mat <- as.matrix(in_r_by_c)

    ## predict cell types
    esti_list <- lapply(seq_len(n_cells), cell_by_type_prob)
    esti_mat <- do.call(rbind(esti_list))
    rownames(esti_mat) <- colnames(in_r_by_c)
    colnames(esti_mat) <- cell_type_names
  }

  return(esti_mat)
}

#' Estimate cell type label for new data using all peaks as input
#'
#' We will identify relevant peaks
#'
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
  ## filter peaks to keep only relevant ones
  r_by_t <- r_by_t[pks_sel, ]
  in_r_by_c <- in_r_by_c[pks_sel, ]

  esti_mat <- estimate_label_default(
    r_by_t = r_by_t, in_r_by_c = in_r_by_c, alpha = alpha
  )

  return(esti_mat)
}
