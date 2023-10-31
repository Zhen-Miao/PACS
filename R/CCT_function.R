### cauchy combination test function


#' Cauchy Combination Test accelerated
#'
#' This is adopted from the Liu and Xie JASA 2020 publication and their GitHub
#' page. Please also cite their initial publication
#' if you used this part of codes. Thaks!
#'
#' @param pval_mat A matrix of p values where each row is a feature and
#'  each column is a p value from one individual test, so that we want to
#'  compute aggregated p values for each row
#'
#' @return A vector of p value for each feature after Cauchy combination
#' @export
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
  pval <- 1 - stats::pcauchy(cct.stat)
  if (sum(is.large_stat) != 0) {
    pval[is.large_stat] <- (1 / cct.stat[is.large_stat]) / pi
  }

  return(pval)
}






#' Cauchy combination test internal (no input checking)
#'
#' This is adopted from the Liu and Xie JASA 2020 publication and their GitHub
#' page. Please also cite their initial publication
#' if you used this part of codes. Thaks!
#'
#' @param pvals A vector of p values to be combined
#'
#' @return A single numeric value of p value after combination
#' @noRd
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
    pval <- 1 - stats::pcauchy(cct.stat)
  }
  return(pval)
}


