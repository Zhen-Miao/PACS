## wrapper for Seurat

## make sure set default null operator exists
if (exists("%||%", envir = baseenv())) {
  `%||%` <- get("%||%", envir = baseenv())
} else {
  `%||%` <- function(x, y) {
    if (is_null(x)) y else x
  }
}

# FindMarkers helper function for cell grouping error checking
ValidateCellGroups_Seurat <- function(
    object,
    cells.1,
    cells.2,
    min.cells.group
) {
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
}





#' Differential test with PACS on Seurat object
#' @description
#' This function performs differential testing on a Seurat object.
#'  Note: Seurat must be installed to use this function.
#' @param cells.1 Vector of cell names belonging to group 1
#' @param cells.2 Vector of cell names belonging to group 2
#' @param features Genes to test. Default is to use all genes
#' @param slot Slot to pull data from; note that if \code{test.use} is
#' "negbinom", "poisson", or "DESeq2", \code{slot} will be set to "counts"
#' @param logfc.threshold Limit testing to genes which show, on average, at least
#' X-fold difference (log-scale) between the two groups of cells. Default is 0.1
#' Increasing logfc.threshold speeds up the function, but can miss weaker signals.
#' If the \code{slot} parameter is "scale.data" no filtering is performed.
#' @param min.pct  only test genes that are detected in a minimum fraction of
#' min.pct cells in either of the two populations. Meant to speed up the function
#' by not testing genes that are very infrequently expressed. Default is 0.01
#' @param min.diff.pct  only test genes that show a minimum difference in the
#' fraction of detection between the two groups. Set to -Inf by default
#' @param verbose Print a progress bar once expression testing begins
#' @param only.pos Only return positive markers (FALSE by default)
#' @param max.cells.per.ident Down sample each identity class to a max number.
#' Default is no downsampling. Not activated by default (set to Inf)
#' @param random.seed Random seed for downsampling
#' @param latent.vars Variables to test
#' @param min.cells.feature Minimum number of cells expressing the feature
#' in at least one of the two groups,
#' currently only used for poisson and negative binomial tests
#' @param min.cells.group Minimum number of cells in one of the groups
#' @param fc.results data.frame from FoldChange
#' @param densify Convert the sparse matrix to a dense form before running the
#' DE test. This can provide speedups but might require higher memory;
#' default is FALSE
#'
#'
#' @importFrom Matrix rowMeans
#' @importFrom stats p.adjust
#'
#' @rdname FindMarkersPACS
#' @concept differential_expression
#' @noRd
#' @method FindMarkersPACS default
#'
.FindMarkersPACS_internal <- function(
    object,
    slot = "data",
    cells.1 = NULL,
    cells.2 = NULL,
    features = NULL,
    logfc.threshold = 0.1,
    min.pct = 0.01,
    min.diff.pct = -Inf,
    verbose = TRUE,
    only.pos = FALSE,
    max.cells.per.ident = Inf,
    random.seed = 1,
    latent.vars = NULL,
    min.cells.feature = 3,
    min.cells.group = 3,
    fc.results = NULL,
    densify = FALSE,
    ...
) {

  ## make sure package is available
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop(paste("Seurat must be installed to use this function.",
         "Please install it using install.packages('Seurat')."))
  }

  ## validate input parameters
  ValidateCellGroups_Seurat(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )

  ## set features
  features <- features %||% rownames(x = object)

  # feature selection (based on percentages)
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(x = alpha.min) <- rownames(x = fc.results)
  features <- names(x = which(x = alpha.min >= min.pct))
  if (length(x = features) == 0) {
    warning("No features pass min.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features <- names(
    x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
  )
  if (length(x = features) == 0) {
    warning("No features pass min.diff.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }

  # feature selection (based on logFC)
  if (slot != "scale.data") {
    total.diff <- fc.results[, 1] #first column is logFC
    names(total.diff) <- rownames(fc.results)
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff >= logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) >= logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      warning("No features pass logfc.threshold threshold; returning empty data.frame")
      return(fc.results[features, ])
    }
  }

  # subsample cell groups if they are too large
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }
  if (inherits(x = object, what = "IterableMatrix")){
    if(test.use != "wilcox"){
      stop("Differential expression with BPCells currently only supports the 'wilcox' method.",
           " Please rerun with test.use = 'wilcox'")
    }
    data.use <- object[features, c(cells.1, cells.2), drop = FALSE]
    groups <- c(rep("foreground", length(cells.1)), rep("background", length(cells.2)))
    de.results <- suppressMessages(
      BPCells::marker_features(data.use, group = groups, method = "wilcoxon")
    )
    de.results <- subset(de.results, de.results$foreground == "foreground")
    de.results <- data.frame(feature = de.results$feature,
                             p_val = de.results$p_val_raw)
    rownames(de.results) <- de.results$feature
    de.results$feature <- NULL
  } else {
    de.results <- PerformDE(
      object = object,
      cells.1 = cells.1,
      cells.2 = cells.2,
      features = features,
      test.use = test.use,
      verbose = verbose,
      min.cells.feature = min.cells.feature,
      latent.vars = latent.vars,
      densify = densify,
      ...
    )
  }
  de.results <- cbind(de.results, fc.results[rownames(x = de.results), , drop = FALSE])
  if (only.pos) {
    de.results <- de.results[de.results[, 2] > 0, , drop = FALSE]
  }
  if (test.use %in% DEmethods_nocorrect()) {
    de.results <- de.results[order(-de.results$power, -de.results[, 1]), ]
  } else {
    de.results <- de.results[order(de.results$p_val, -abs(de.results$pct.1-de.results$pct.2)), ]
    de.results$p_val_adj = p.adjust(
      p = de.results$p_val,
      method = "bonferroni",
      n = nrow(x = object)
    )
  }
  return(de.results)
}

#' @param fc.slot Slot used to calculate fold-change - will also affect the
#' default for \code{mean.fxn}, see below for more details.
#' @param pseudocount.use Pseudocount to add to averaged expression values when
#' calculating logFC. 1 by default.
#' @param norm.method Normalization method for fold change calculation when
#' \code{slot} is \dQuote{\code{data}}
#' @param mean.fxn Function to use for fold change or average difference calculation.
#' The default depends on the the value of \code{fc.slot}:
#' \itemize{
#'  \item{"counts"} : difference in the log of the mean counts, with pseudocount.
#'  \item{"data"} : difference in the log of the average exponentiated data, with pseudocount.
#'  This adjusts for differences in sequencing depth between cells, and assumes that "data"
#'  has been log-normalized.
#'  \item{"scale.data"} : difference in the means of scale.data.
#' }
#' @param fc.name Name of the fold change, average difference, or custom function column
#' in the output data.frame. If NULL, the fold change column will be named
#' according to the logarithm base (eg, "avg_log2FC"), or if using the scale.data
#' slot "avg_diff".
#' @param base The base with respect to which logarithms are computed.
#'
#' @rdname FindMarkers
#' @concept differential_expression
#' @export
#' @method FindMarkers Assay
#'
FindMarkersPACS <- function(
    object,
    slot = "data",
    cells.1 = NULL,
    cells.2 = NULL,
    features = NULL,
    fc.slot = "data",
    pseudocount.use = 1,
    norm.method = NULL,
    mean.fxn = NULL,
    fc.name = NULL,
    base = 2,
    ...
) {

  if (length(x = Layers(object = object, search = slot)) > 1) {
    stop(slot, " layers are not joined. Please run JoinLayers")
  }
  data.use <-  Seurat::GetAssayData(object = object, slot = data.slot)
  fc.results <- Seurat::FoldChange(
    object = object,
    slot = fc.slot,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    pseudocount.use = pseudocount.use,
    mean.fxn = mean.fxn,
    fc.name = fc.name,
    base = base,
    norm.method = norm.method
  )
  de.results <- .FindMarkersPACS_internal(
    object = data.use,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    test.use = test.use,
    fc.results = fc.results,
    ...
  )
  return(de.results)
}



pacs_test_seurat <- function(object, cell_types_of_interest,
                             by_identity = TRUE,
                             meta_to_keep, formula_full,
                             formula_null, pic_matrix,
                             n_peaks_per_round = NULL,
                             T_proportion_cutoff = 0.2,
                             cap_rates, par_initial_null = NULL,
                             par_initial_full = NULL, n_cores = 1,
                             verbose = TRUE) {

  ## need Seurat package installed
  require("Seurat")

  ## check the object
  if(!is(object, 'Seurat')){
    stop("The input object must be a Seurat object!")
  }

  ## check formula
  obmeta = object_sub@meta.data
  vars_in_formula_full <- all.vars(formula_full)
  vars_in_formula_null <- all.vars(formula_null)
  if(!all(vars_in_formula_full %in% colnames(obmeta)) |
     !all(vars_in_formula_null %in% colnames(obmeta))){
    stop(paste("Not all variables in formula are available in meta.data slot,",
               "Please make sure to include all variables in the meta.data",
               "slot. "))
  }


  object_sub <- subset(object, idents = cell_types_of_interest)



  ## convert everthing into a character -- drop all
  # pbmeta <- droplevels(pbmeta)
  pbmeta[] <- lapply(pbmeta, function(x) if(is.factor(x)) as.character(x) else x)





  return(object)
}
