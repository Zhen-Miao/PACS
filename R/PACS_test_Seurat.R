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
.FindMarkersPACS_default <- function(
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

  ## run differential test
    de.results <- PerformDEPACS(
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

  de.results <- cbind(de.results, fc.results[rownames(x = de.results), , drop = FALSE])
  if (only.pos) {
    de.results <- de.results[de.results[, 2] > 0, , drop = FALSE]
  }

  de.results <- de.results[order(de.results$p_val, -abs(de.results$pct.1-de.results$pct.2)), ]
  de.results$p_val_adj = p.adjust(
    p = de.results$p_val,
    method = "bonferroni",
    n = nrow(x = object)
  )

  return(de.results)
}

#' Differential test with PACS on Seurat object
#' @description
#' Run differential test using the similar style as in Seurat
#' `Seurat::FindMarkers()` function.
#' Note 1: This function is modified from Seurat package v5.1.0
#' Note 2: Seurat must be installed to use this function.
#' @param slot Slot to pull data from; note that for PACS,
#'  it should be set to "counts"
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
#' @rdname FindMarkersPACS
#' @concept differential_expression
#' @export
#' @method FindMarkersPACS Assay
#'
.FindMarkersPACS_assay <- function(
    object,
    slot = "counts",
    cells.1 = NULL,
    cells.2 = NULL,
    features = NULL,
    fc.slot = "counts",
    pseudocount.use = 1,
    norm.method = NULL,
    mean.fxn = NULL,
    fc.name = NULL,
    base = 2,
    ...
) {

  ## for scATAC-seq data, we should use count slot
  data.slot <- "counts"

  if (length(x = Layers(object = object, search = slot)) > 1) {
    stop(slot, " layers are not joined. Please run JoinLayers")
  }

  ## obtain the data matrix
  data.use <-  Seurat::GetAssayData(object = object, slot = data.slot)

  ## calculate FC
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

  ## differential test
  de.results <- .FindMarkersPACS_default(
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


#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
#' @param group.by Regroup cells into a different identity class prior to
#' performing differential expression (see example)
#' @param subset.ident Subset a particular identity class prior to regrouping.
#' Only relevant if group.by is set (see example)
#' @param assay Assay to use in differential expression testing
#' @param reduction Reduction to use in differential expression testing - will
#' test for DE on cell embeddings
#'
#' @rdname FindMarkersPACS
#' @concept differential_expression
#' @export
#' @method FindMarkersPACS Seurat
#'
FindMarkersPACS <- function(
    object,
    ident.1 = NULL,
    ident.2 = NULL,
    latent.vars = NULL,
    group.by = NULL,
    subset.ident = NULL,
    assay = NULL,
    reduction = NULL,
    ...
) {
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  if (length(x = ident.1) == 0) {
    stop("At least 1 ident must be specified in `ident.1`")
  }

  # select which data to use
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(x = data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(x = data.use)
  }

  cells <- IdentsToCells(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    cellnames.use = cellnames.use
  )
  cells <- sapply(
    X = cells,
    FUN = intersect,
    y = cellnames.use,
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  if (!all(vapply(X = cells, FUN = length, FUN.VALUE = integer(length = 1L)))) {
    abort(
      message = "Cells in one or both identity groups are not present in the data requested"
    )
  }

  # fetch latent.vars
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(cells$cells.1, cells$cells.2)
    )
  }

  # check normalization method
  norm.command <- paste0("NormalizeData.", assay)
  norm.method <- if (norm.command %in% Command(object = object) && is.null(x = reduction)) {
    Command(
      object = object,
      command = norm.command,
      value = "normalization.method"
    )
  } else if (length(x = intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"), y = Command(object = object)))) {
    command <- intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"), y = Command(object = object))[1]
    Command(
      object = object,
      command = command,
      value = "normalization.method"
    )
  } else {
    NULL
  }

  de.results <- .FindMarkersPACS_assay(
    object = data.use,
    latent.vars = latent.vars,
    cells.1 = cells$cells.1,
    cells.2 = cells$cells.2,
    norm.method = norm.method,
    ...
  )

  return(de.results)
}


PerformDEPACS <- function(
    object,
    cells.1,
    cells.2,
    features,
    verbose,
    min.cells.feature,
    latent.vars,
    densify,
    ...
) {

  data.use <- object[features, c(cells.1, cells.2), drop = FALSE]
  if (densify){
    data.use <- as.matrix(x = data.use)
  }
  de.results <-
    LRDETest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    )
  return(de.results)
}

