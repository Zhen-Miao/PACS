## wrapper for covariate matrix


#'  Profile likelihood ratio test with PACS
#'
#' @importFrom stats model.matrix
#'
#' @param covariate_meta.data A data.frame with columns representing the
#'   covariates and rows representing cells
#' @param formula_full A formula object representing the full model. For
#'   example, ~ cell_type + batch
#' @param formula_null A formula object representing the null model. For
#'   example, ~ batch
#' @param pic_matrix The input region-by-cell PIC matrix
#' @param cap_rates A vector of capturing probability for each cell
#' @param par_initial_null Initialized values of estimated parameters for the
#'   null model, we do not
#'   need to specify unless there are reasons to do so. Default = NULL
#' @param par_initial_full Initialized values of estimated parameters for the
#'   null model, we do not
#'   need to specify unless there are reasons to do so. Default = NULL
#' @param n_cores number of cores for multi-core computation
#'
#' @return A list of two elements, pacs_converged is a vector
#'   of length 2*n_peaks
#'   representing the convergence status of the peak in the null and full model
#'   and pacs_p_val is a vector of length n_peaks representing the p values for
#'   each peak.
#' @export
#'
pacs_test_logit <- function(covariate_meta.data, formula_full,
                            formula_null, pic_matrix,
                            cap_rates, par_initial_null = NULL,
                            par_initial_full = NULL, n_cores = 1) {
  ### construct model matrix
  X_full <- model.matrix(formula_full, data = covariate_meta.data)
  X_null <- model.matrix(formula_null, data = covariate_meta.data)

  ### re-order the covariate matrix so that the last few columns are
  ### parameters of interest
  pars_of_interest <- setdiff(colnames(X_full), colnames(X_null))
  X_full <- X_full[, c(colnames(X_null), pars_of_interest)]
  index_poi <- which(colnames(X_full) %in% pars_of_interest)

  ### number of parameters in both matrices
  n_para_full <- dim(X_full)[2]

  ### initialize the estimated parameters
  if (is.null(par_initial_full)) {
    par_initial_full <- rep.int(0.05, n_para_full)
  }else if (length(par_initial_full) != n_para_full) {
    par_initial_full <- rep(par_initial_full[1], times = n_para_full)
  }

  if (is.null(par_initial_null)) {
    par_initial_null <- par_initial_full
  }else if (length(par_initial_null) != n_para_full) {
    par_initial_null <- rep(par_initial_null[1], times = n_para_full)
  }

  ## make sure the null part is set to zero
  par_initial_null[!(colnames(X_full) %in% colnames(X_null))] <- 0

  ### fit the null model
  null_para <- estimate_parameters_null(
    r_by_c = pic_matrix,
    design_mat = X_full, ## full design mat
    par_initial = par_initial_null,
    hold_zero = index_poi,
    cap_rate_vec = cap_rates,
    mc.cores = n_cores
  )

  ### fit the full model
  full_para <- estimate_parameters(
    r_by_c = pic_matrix,
    design_mat = X_full,
    par_initial = par_initial_full,
    cap_rate_vec = cap_rates,
    mc.cores = n_cores
  )

  ### converge check and prepare for test
  conv_mat <- c(
    null_para[dim(null_para)[1], ],
    full_para[dim(full_para)[1], ]
  ) ## convergence vector
  null_para <- null_para[1:(dim(null_para)[1] - 1),, drop = FALSE ]
  full_para <- full_para[1:(dim(full_para)[1] - 1),, drop = FALSE ]

  ### p value
  pacs_p_val <- compare_models(
    x_full = X_full,
    theta_estimated_full = full_para,
    x_null = X_full,
    theta_estimated_null = null_para,
    q_vec = cap_rates,
    c_by_r = t(pic_matrix),
    df_test = length(index_poi),
    mc.cores = n_cores
  )

  names(pacs_p_val) <- rownames(pic_matrix)


  return(list(pacs_converged = conv_mat, pacs_p_val = pacs_p_val))
}
