### cauchy combination test function


CCT_internal_horizontal <- function(pval_mat) {
  ## each row is a feature, each column is one p value to be integrated

  if (is.null(rownames(pval_mat))) {
    rownames(pval_mat) <- paste("f", seq_len(nrow(pval_mat)), sep = "_")
  }

  if (is.null(colnames(pval_mat))) {
    colnames(pval_mat) <- paste("p", seq_len(ncol(pval_mat)), sep = "_")
  }

  rname_ori <- rownames(pval_mat)



  weights <- matrix(1 / ncol(pval_mat), ncol = ncol(pval_mat),
                    nrow = nrow(pval_mat))

  is.small <- rowSums(pval_mat < 1e-16) > 0

  if (sum(is.small) != 0) {
    ## non-small values
    pval_n_small <- pval_mat[!is.small, ]
    weights_n_small <- weights[!is.small, ]
    cct.stat_n_small <- rowSums(weights_n_small * tan((0.5 - pval_n_small) * pi))
    names(cct.stat_n_small) <- rownames(pval_n_small)

    ## small values
    pval_small <- pval_mat[is.small, ]

    which.small <- pval_small < 1e-16
    cct.stat_small <- vector(length = nrow(pval_small))

    for (i in seq_len(nrow(pval_small))) {
      w_small <- which.small[i, ]
      cct.stat_small[i] <-
        sum((weights[i, w_small] / pval_small[i, w_small]) / pi)
      cct.stat_small[i] <- cct.stat_small[i] +
        sum(weights[i, !w_small] * tan((0.5 - pval_small[i, !w_small]) * pi))
    }
    names(cct.stat_small) <- rownames(pval_small)

    ## combine small and large
    cct.stat <- c(cct.stat_n_small, cct.stat_small)
    cct.stat <- cct.stat[rname_ori]
  } else {
    pval_n_small <- pval_mat
    cct.stat <- rowSums(weights * tan((0.5 - pval_n_small) * pi))
  }

  is.large_stat <- cct.stat > 1e+15
  pval <- 1 - pcauchy(cct.stat)
  if (sum(is.large_stat) != 0) {
    pval[is.large_stat] <- (1 / cct.stat[is.large_stat]) / pi
  }

  return(pval)
}






CCT_internal <- function(pvals) {
  weights <- rep(1 / length(pvals), length(pvals))

  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
  } else {
    cct.stat <- sum((weights[is.small] / pvals[is.small]) / pi)
    cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
  }

  #### check if the test statistic is very large.
  if (cct.stat > 1e+15) {
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - pcauchy(cct.stat)
  }
  return(pval)
}


CCT <- function(pvals, weights = NULL) {
  #### check if there is NA
  if (sum(is.na(pvals)) > 0) {
    stop("Cannot have NAs in the p-values!")
  }

  #### check if all p-values are between 0 and 1
  if ((sum(pvals < 0) + sum(pvals > 1)) > 0) {
    stop("All p-values must be between 0 and 1!")
  }

  #### check if there are p-values that are either exactly 0 or 1.
  is.zero <- (sum(pvals == 0) >= 1)
  is.one <- (sum(pvals == 1) >= 1)
  if (is.zero && is.one) {
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero) {
    return(0)
  }
  if (is.one) {
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### check the validity of weights (default: equal weights)
  #### and standardize them.
  if (is.null(weights)) {
    weights <- rep(1 / length(pvals), length(pvals))
  } else if (length(weights) != length(pvals)) {
    stop("The length of weights should be the same as that of the p-values!")
  } else if (sum(weights < 0) > 0) {
    stop("All the weights must be positive!")
  } else {
    weights <- weights / sum(weights)
  }

  #### check if there are very small non-zero p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0) {
    cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
  } else {
    cct.stat <- sum((weights[is.small] / pvals[is.small]) / pi)
    cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
  }

  #### check if the test statistic is very large.
  if (cct.stat > 1e+15) {
    pval <- (1 / cct.stat) / pi
  } else {
    pval <- 1 - pcauchy(cct.stat)
  }
  return(pval)
}
