## wrapper for covariate matrix


#' likelihood ratio test with PACS -- cumulative logit
#'
#' @importFrom stats model.matrix
#' @param covariate_meta.data A data.frame with columns representing the
#'   covariates and rows representing cells
#' @param formula_full A formula object representing the full model. For
#'   example, ~ cell_type + batch
#' @param formula_null A formula object representing the null model. For
#'   example, ~ batch
#' @param pic_matrix The input region-by-cell PIC matrix
#' @param max_T The maximum accessibility category considered, default = 2.
#'   For `method = "exact"`, observed counts greater than `max_T` are
#'   top-coded. Under `capture = "B"` this is exact and the final category
#'   means `max_T` or more. Under `capture = "A"` it is an approximation:
#'   capping and fragment thinning do not commute, so `max_T` acts as a cap on
#'   the *latent* count and should be set high enough that observed counts
#'   above it are rare. A warning reports how many were top-coded.
#' @param cap_rates A vector of capturing probability for each cell
#' @param par_initial_null Initialized values of estimated parameters for the
#'   null model, we do not
#'   need to specify unless there are reasons to do so. Default = NULL
#' @param par_initial_full Initialized values of estimated parameters for the
#'   null model, we do not
#'   need to specify unless there are reasons to do so. Default = NULL
#' @param n_cores number of cores for multi-core computation
#' @param method One of `"stacked"` (default, back-compat) or `"exact"`.
#'   `"stacked"` uses the original stack-and-treat-as-binary approximation.
#'   `"exact"` uses the proper cumulative-logit likelihood with capture-rate
#'   correction (`notes/cumulative_logit_math.md`) and an unpenalized
#'   maximum-likelihood fit. Exact-path p-values are ordinary
#'   likelihood-ratio tests and are `NA` for non-converged or boundary fits.
#' @param capture Capture model for `method = "exact"`; ignored by
#'   `method = "stacked"`. `"B"` (default, back-compat) is cell-level
#'   all-or-nothing dropout: with probability `q_i` the cell is observed
#'   perfectly, otherwise its count collapses to zero. `"A"` is fragment-level
#'   thinning, `M_i | Y_i = k ~ Binomial(k, q_i)`, which is the natural
#'   generalisation of the binary PACS capture model and the appropriate
#'   choice when a partially captured cell should still be able to yield a
#'   nonzero count. The two coincide at `max_T = 1` and diverge as counts of
#'   2 or more become common. See Options A and B in
#'   `notes/cumulative_logit_math.md`.
#' @details The unpenalized exact MLE can fail to exist or have singular
#'   information for very sparse peaks. In review simulations, the exact-path
#'   `NA` rate rose from 0% for dense peaks to 3% for moderately sparse peaks
#'   (`n = 300`, `alpha = c(-2.5, -4)`), 38% for still sparser peaks
#'   (`n = 300`, `alpha = c(-3.5, -5)`), and 46% with fewer cells
#'   (`n = 100`, `alpha = c(-2.5, -4)`). These rates are scenario-specific,
#'   not general guarantees. Withholding a p-value is a conservative failure
#'   policy for the affected peak—it avoids turning a failed fit into a false
#'   positive—but the resulting loss of analyzable peaks reduces power.
#'
#' @return A list of two elements. `pacs_converged` has length
#'   `2 * n_peaks`, with null-fit statuses followed by full-fit statuses. For
#'   the exact path, status 1 means converged, 2 means a singular or non-finite
#'   scoring system, 3 means the iteration limit was reached, 4 means
#'   step-halving found no acceptable update, and 5 means the supplied starting
#'   value had a non-finite log-likelihood. `pacs_p_val` contains one p-value per
#'   peak; exact-path inference is `NA` unless both statuses are 1 and both fits
#'   are interior.
#' @export
#'
pacs_test_cumu <- function(covariate_meta.data, formula_full,
                           formula_null, pic_matrix, max_T = 2,
                           cap_rates, par_initial_null = NULL,
                           par_initial_full = NULL, n_cores = 1,
                           method = c("stacked", "exact"),
                           capture = c("B", "A")) {
  method <- match.arg(method)
  capture_supplied <- !missing(capture)
  capture <- match.arg(capture)
  if (method != "exact" && capture_supplied) {
    ## Silently ignoring a capture-model choice would misrepresent what was
    ## fitted, so say so rather than let it pass unnoticed.
    warning(
      "`capture` only applies to method = \"exact\" and was ignored.",
      call. = FALSE
    )
  }
  if (method == "exact") {
    return(pacs_test_cumu_exact(
      covariate_meta.data = covariate_meta.data,
      formula_full = formula_full,
      formula_null = formula_null,
      pic_matrix = pic_matrix,
      max_T = max_T,
      cap_rates = cap_rates,
      par_initial_null = par_initial_null,
      par_initial_full = par_initial_full,
      n_cores = n_cores,
      capture = capture
    ))
  }

  ### construct model matrix
  X_full <- model.matrix(formula_full, data = covariate_meta.data)
  X_null <- model.matrix(formula_null, data = covariate_meta.data)

  ### stack the model matrix
  X_full_beta <- X_full[, 2:ncol(X_full), drop = FALSE]
  X_full_beta <- replicate(max_T, X_full_beta, simplify = FALSE)
  X_full_stacked <- do.call(rbind, X_full_beta)
  colnames(X_full_stacked) <- colnames(X_full)[2:ncol(X_full)]

  ## the intersection part
  A <- diag(nrow = max_T, ncol = max_T)
  X_alpha <- A[rep(seq_len(nrow(A)), each = nrow(X_full)), ]
  colnames(X_alpha) <- paste("intersection_", 1:max_T, sep = "")

  ## construct the stacked X_full
  X_full <- cbind(X_alpha, X_full_stacked)

  if (ncol(X_null) != 1) {
    X_null_beta <- X_null[, 2:ncol(X_null), drop = FALSE]
    X_null_beta <- replicate(max_T, X_null_beta, simplify = FALSE)
    X_null_stacked <- do.call(rbind, X_null_beta)
    colnames(X_null_stacked) <- colnames(X_null)[2:ncol(X_null)]

    ## construct the stacked X_null
    X_null <- cbind(X_alpha, X_null_stacked)
  } else {
    X_null <- X_alpha
  }




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

  ## also stack the input matrix, and set the values accordingly
  pic_cumu <- rep(list(), length = max_T)
  for (t in 1:max_T) {
    pic_cumu[[t]] <- ifelse(pic_matrix >= t, 1, 0)
  }
  pic_matrix <- do.call(cbind, pic_cumu)

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
  null_para <- null_para[1:(dim(null_para)[1] - 1),, drop = FALSE  ]
  full_para <- full_para[1:(dim(full_para)[1] - 1),, drop = FALSE  ]

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


## Exact-likelihood path. Internal driver invoked from pacs_test_cumu when
## method = "exact". See notes/cumulative_logit_math.md for the model;
## fitting routines live in R/param_estimate_cumu.R.
pacs_test_cumu_exact <- function(covariate_meta.data, formula_full,
                                 formula_null, pic_matrix, max_T,
                                 cap_rates, par_initial_null,
                                 par_initial_full, n_cores,
                                 capture = c("B", "A")) {
  capture <- match.arg(capture)
  if (length(max_T) != 1L || !is.finite(max_T) || max_T < 1L ||
      max_T != as.integer(max_T)) {
    stop("max_T must be a positive integer.", call. = FALSE)
  }
  T <- as.integer(max_T)
  M_dense <- as.matrix(pic_matrix)
  if (any(!is.finite(M_dense)) || any(M_dense < 0) ||
      any(M_dense != floor(M_dense))) {
    stop(
      "pic_matrix must contain finite, non-negative integer counts.",
      call. = FALSE
    )
  }
  if (nrow(covariate_meta.data) != ncol(M_dense)) {
    stop(
      "covariate_meta.data must have one row per column of pic_matrix.",
      call. = FALSE
    )
  }
  if (length(cap_rates) != ncol(M_dense) ||
      any(!is.finite(cap_rates)) || any(cap_rates <= 0) ||
      any(cap_rates > 1)) {
    stop(
      "cap_rates must contain one finite value in (0, 1] per cell.",
      call. = FALSE
    )
  }

  ## Under Option B, top-coding the observed count is exact: M_i = Y_i * B_i
  ## with B_i Bernoulli, so min(M_i, T) = min(Y_i, T) * B_i is again an Option
  ## B observation with latent min(Y_i, T). Category T legitimately means
  ## "Y_i >= T".
  ##
  ## Under Option A it is not. Capping and thinning do not commute:
  ## min(Binom(Y, q), T) is not Binom(min(Y, T), q) -- at Y = 5, q = 0.5,
  ## T = 2 the two give Pr(M = 2) = 0.81 and 0.25. So under Option A the cap
  ## on the latent Y is a modelling assumption ("no cell carries more than T
  ## fragments"), not a relabelling of the top category, and max_T must be
  ## large enough that exceedances are rare.
  exceedances <- sum(M_dense > T)
  if (capture == "A" && exceedances > 0L) {
    warning(sprintf(
      paste0(
        "%d observed count(s) exceed max_T = %d and were top-coded. Under ",
        "capture = \"A\" the latent count is capped at max_T, and capping ",
        "does not commute with fragment thinning, so this is an ",
        "approximation rather than a relabelling of the top category. ",
        "Consider raising max_T."
      ),
      exceedances, T
    ), call. = FALSE)
  }
  pic_matrix <- pmin(M_dense, T)
  dimnames(pic_matrix) <- dimnames(M_dense)

  X_full <- model.matrix(formula_full, data = covariate_meta.data)
  X_null <- model.matrix(formula_null, data = covariate_meta.data)

  ## Strip leading intercept columns (alpha plays the intercept role per
  ## threshold). Both formulas should produce a leading "(Intercept)" column.
  if ("(Intercept)" %in% colnames(X_full)) {
    X_full <- X_full[, setdiff(colnames(X_full), "(Intercept)"), drop = FALSE]
  }
  if ("(Intercept)" %in% colnames(X_null)) {
    X_null <- X_null[, setdiff(colnames(X_null), "(Intercept)"), drop = FALSE]
  }

  ## Identify parameters of interest (in full but not null) and re-order so
  ## they sit at the end of the design matrix.
  pars_of_interest <- setdiff(colnames(X_full), colnames(X_null))
  if (length(pars_of_interest) == 0L) {
    stop("formula_full and formula_null have the same covariates; nothing to test.")
  }
  X_full <- X_full[, c(colnames(X_null), pars_of_interest), drop = FALSE]
  p_beta <- ncol(X_full)

  ## hold_zero indices in the full theta = (atilde_1..T, beta_1..p):
  ## thresholds always free; betas of interest set to zero under null.
  beta_poi_idx <- which(colnames(X_full) %in% pars_of_interest)
  hold_zero <- T + beta_poi_idx

  ## Initial values: warm-start atilde from the entry-pooled empirical
  ## marginals of pic_matrix; beta = 0. Single shared starting point for
  ## all peaks (matches the binary path's `par_initial = rep.int(0.05, n)`
  ## convention).
  if (is.null(par_initial_full) || length(par_initial_full) != T + p_beta) {
    M_pool <- as.numeric(pic_matrix)
    par_initial_full <- warm_start_theta(
      M = M_pool, q = cap_rates, T = T, p_beta = p_beta
    )
  }
  if (is.null(par_initial_null)) {
    par_initial_null <- par_initial_full
  } else if (length(par_initial_null) != T + p_beta) {
    par_initial_null <- par_initial_full
  }
  par_initial_null[hold_zero] <- 0

  null_para <- estimate_parameters_cumu_null(
    r_by_c = pic_matrix, X = X_full, theta_initial = par_initial_null,
    hold_zero = hold_zero, q_vec = cap_rates, T = T, mc.cores = n_cores,
    capture = capture
  )
  full_para <- estimate_parameters_cumu(
    r_by_c = pic_matrix, X = X_full, theta_initial = par_initial_full,
    q_vec = cap_rates, T = T, mc.cores = n_cores, capture = capture
  )

  conv_null <- null_para[nrow(null_para), ]
  conv_full <- full_para[nrow(full_para), ]
  conv_mat <- c(conv_null, conv_full)
  null_para <- null_para[seq_len(nrow(null_para) - 1L), , drop = FALSE]
  full_para <- full_para[seq_len(nrow(full_para) - 1L), , drop = FALSE]

  c_by_r <- t(pic_matrix)

  pacs_p_val <- compare_models_cumu(
    x_full = X_full,
    theta_estimated_full = full_para,
    theta_estimated_null = null_para,
    q_vec = cap_rates, c_by_r = c_by_r, T = T,
    df_test = length(beta_poi_idx),
    conv_full = conv_full, conv_null = conv_null,
    capture = capture, mc.cores = n_cores
  )
  names(pacs_p_val) <- rownames(pic_matrix)

  list(pacs_converged = conv_mat, pacs_p_val = pacs_p_val)
}
