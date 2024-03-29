## param_estimate_simple

##################################################################
## The functions in this file are only for method comparisons
## we did not use functions here in the PACS workflow
##################################################################

#' Add a pseudo-count to data matrix
#'
#' Please note: This function is only for method comparison purposes,
#' we did not use this function in the PACS workflow
#'
#'
#' @param r_by_c region by cell matrix
#' @param group_set A set of cell type labels
#' @param group_labels A vector of cell type labels for each cell
#' @param batch_labels A vector of batch labels for each cell
#' @param capturing_rate A vector of capturing rates for each cell
#'
#' @return A list with three elements, r_by_c_pseudo, group_labels_pseudo,
#'   and capturing_rate_pseudo
#'
#' @noRd
add_pseudo_count <- function(
    r_by_c,
    group_set,
    group_labels,
    batch_labels,
    capturing_rate) {
  if (!is.na(batch_labels)) {
    n_batch <- length(unique(batch_labels))
    batch_set <- unique(batch_labels)
  } else {
    n_batch <- 1
  }

  ## get input info
  n_group <- length(group_set)
  # n_peak <- dim(r_by_c)[1]
  # n_cell <- dim(r_by_c)[2]

  if (!is.na(batch_labels)) {
    new_batch_labels <- rep(batch_set, times = n_group, each = 2)
  }


  ## create new batch and group labels
  new_group_labels <- rep(group_set, times = 1, each = (2 * n_batch))
  new_capturing_rate <- rep(0.95, times = (2 * n_batch * n_group))

  ## create new matrix
  new_r_by_c <- matrix(
    data = 0, nrow = nrow(r_by_c),
    ncol = (2 * n_batch * n_group)
  )
  for (i in 1:(n_batch * n_group)) {
    new_r_by_c[, 2 * i] <- 1
  }

  ## append the values
  r_by_c_pseudo <- cbind(r_by_c, new_r_by_c)

  group_labels_pseudo <- append(group_labels, new_group_labels)
  capturing_rate_pseudo <- append(capturing_rate, new_capturing_rate)
  if (!is.na(batch_labels)) {
    batch_labels_pseudo <- append(batch_labels, new_batch_labels)
  }

  ## return values
  if (!is.na(batch_labels)) {
    ret <- list(
      r_by_c_pseudo = r_by_c_pseudo,
      batch_labels_pseudo = batch_labels_pseudo,
      group_labels_pseudo = group_labels_pseudo,
      capturing_rate_pseudo = capturing_rate_pseudo
    )
  } else {
    ret <- list(
      r_by_c_pseudo = r_by_c_pseudo,
      group_labels_pseudo = group_labels_pseudo,
      capturing_rate_pseudo = capturing_rate_pseudo
    )
  }
  return(ret)
}


#' get_r_by_ct_mat
#'
#' Please note: This function is only for method comparison purposes,
#' we did not use this function in the PACS workflow
#'
#' @param r_by_c region by cell matrix
#' @param group_set A set of cell type labels
#' @param group_labels A vector of cell type labels for each cell
#' @param batch_labels A vector of batch labels for each cell
#' @param capturing_rate A vector of capturing rates for each cell
#' @param max_p Maximum value for the open probability, default = 0.999
#' @param min_p Minimum value for the open probability, default = 0.0001
#' @param adding_doublet Whether to consider doublet
#' @param doublet_model The model of doublet
#' @param adding_pseudo_count Whether to add pseudo counts
#'
#' @return summary of region by cell type matrix
#'
#' @noRd
get_r_by_ct_mat <- function(
    r_by_c,
    group_set,
    group_labels,
    capturing_rate,
    max_p = 0.999,
    min_p = 0.0001,
    adding_doublet = FALSE,
    doublet_model = "max",
    adding_pseudo_count = TRUE) {
  ## check input peak by cell
  if (dim(r_by_c)[2] != length(group_labels)) {
    stop("n_cells in pbyc should equals to n_cells in labels")
  } else if (dim(r_by_c)[2] != length(capturing_rate)) {
    stop("n_cells in pbyc should equals to n_cells in capturing rate")
  }

  ## check input doublet_model
  if (length(doublet_model) != 1) {
    stop("please specify only one doublet model at a time")
  } else if (!(doublet_model %in% c("expected", "half", "max", "ave"))) {
    stop("please specify one of the models: expected, half, max, ave")
  }


  pbyt <- matrix(nrow = dim(r_by_c)[1], ncol = length(group_set))
  rownames(pbyt) <- rownames(r_by_c)
  colnames(pbyt) <- group_set

  ## add pesudo_count
  if (adding_pseudo_count) {
    ret_list <- add_pseudo_count(
      r_by_c = r_by_c, group_set = group_set, group_labels = group_labels,
      batch_labels = NA, capturing_rate = capturing_rate
    )

    r_by_c <- ret_list$r_by_c_pseudo
    group_labels <- ret_list$group_labels_pseudo
    capturing_rate <- ret_list$capturing_rate_pseudo
  }

  sum_capturing_rate <- vector(length = length(group_set))
  names(sum_capturing_rate) <- group_set
  for (i in group_set) {
    sum_capturing_rate[i] <- sum(capturing_rate[group_labels == i])
  }
  for (i in group_set) {
    # print(paste(i,'running'))
    pbyt_vec <- (Matrix::rowSums(r_by_c[, group_labels == i])) / sum_capturing_rate[i]

    ## correct for values >=1 and =0
    pbyt_vec[pbyt_vec > max_p] <- max_p
    pbyt_vec[pbyt_vec < min_p] <- min_p
    pbyt[, i] <- pbyt_vec
  }
  if (adding_doublet == TRUE) {
    ### now we do all doublets together

    for (i in seq_along(group_set)) {
      for (j in seq_along(group_set)) {
        ## require different cell types
        if (i >= j) {
          next
        } else {
          cell_type_sample_i <- group_set[i]
          cell_type_sample_j <- group_set[j]
        }

        ## combine cell types
        if (doublet_model == "expected") {
          doublet <- 1 - ((1 - pbyt[, cell_type_sample_i]) *
            (1 - pbyt[, cell_type_sample_j]))
        } else if (doublet_model == "half") { ## not recommended
          doublet <- 0.5 * (1 - ((1 - pbyt[, cell_type_sample_i]) *
            (1 - pbyt[, cell_type_sample_j])))
        } else if (doublet_model == "max") { ## not recommended
          doublet <- pmax(
            pbyt[, cell_type_sample_i],
            pbyt[, cell_type_sample_j]
          )
        } else if (doublet_model == "ave") { ### not recommended
          doublet <- (0.5 * pbyt[, cell_type_sample_i]) +
            (0.5 * pbyt[, cell_type_sample_j])
        }

        doublet[doublet > max_p] <- max_p
        doublet[doublet < min_p] <- min_p
        pbyt <- cbind(pbyt, doublet)
        colnames(pbyt)[dim(pbyt)[2]] <- paste("doublet_", cell_type_sample_i,
          cell_type_sample_j,
          sep = ""
        )
      }
    }
  }
  return(pbyt)
}


#' Estimate cell type label with capturing rate
#'
#' Please note: This function is only for method comparison purposes,
#' we did not use this function in the PACS workflow
#'
#' @import Matrix
#' @param r_by_t region by type matrix
#' @param in_r_by_c input region by cell matrix
#' @param capturing_rate Capturing probability (rate)
#' @param alpha alpha value for considering negative peaks
#'
#' @return A matrix of cell types by cell matrix, with element representing
#'  probability of being in that cell type
#' @noRd
estimate_label_w_capturing_rate <- function(
    r_by_t,
    in_r_by_c,
    capturing_rate,
    alpha = 1) {
  ## build prediction matrix
  esti_m5 <- matrix(nrow = dim(in_r_by_c)[2], ncol = dim(r_by_t)[2])
  ## row cell, column cell types
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

  if (as.numeric(dim(in_r_by_c)[1]) * as.numeric(dim(in_r_by_c)[2]) >
    (2^31 - 1)) {
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
        pqbyt <- r_by_t * capturing_rate[i_cell]
        lg_pqbyt <- log(pqbyt)
        lg_pqbyt_q <- log1p(-1 * pqbyt)
        # unparallelized
        esti_m5[i_cell + (i - 1) * dim2, ] <-
          colSums(lg_pqbyt * xvec + alpha * lg_pqbyt_q * (1 - xvec))
      }
    }
  } else {
    in_r_by_c_mat <- as.matrix(in_r_by_c)

    ## predict cell types

    # unparallelized
    for (i_cell in seq_len(dim(esti_m5)[1])) {
      pqbyt <- r_by_t * capturing_rate[i_cell]
      lg_pqbyt <- log(pqbyt)
      lg_pqbyt_q <- log1p(-1 * pqbyt)
      # unparallelized
      xvec <- in_r_by_c_mat[, i_cell]
      esti_m5[i_cell, ] <- colSums(lg_pqbyt * xvec + lg_pqbyt_q * (1 - xvec))
    }
    return(esti_m5)
  }
}
