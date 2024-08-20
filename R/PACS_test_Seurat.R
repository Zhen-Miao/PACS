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
.ValidateCellGroups_Seurat <- function(
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
#' @param slot Slot to pull data from
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
#' @importFrom Matrix rowMeans
#' @importFrom stats p.adjust
#'
#' @concept differential_expression
#' @noRd
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
  .ValidateCellGroups_Seurat(
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
#' @concept differential_expression
#' @noRd
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

  if (length(x = SeuratObject::Layers(object = object, search = slot)) > 1) {
    stop(slot, " layers are not joined. Please run JoinLayers")
  }

  ## obtain the data matrix
  data.use <-  SeuratObject::GetAssayData(object = object, slot = data.slot)

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
    fc.results = fc.results,
    ...
  )
  return(de.results)
}

# Function to get all the descendants on a tree of a given node
#
# @param tree Tree object (from ape package)
# @param node Internal node in the tree
#
# @return Returns all descendants of the given node
#
.GetDescendants <- function(tree, node, curr = NULL) {
  if (is.null(x = curr)) {
    curr <- vector()
  }
  daughters <- tree$edge[which(x = tree$edge[, 1] == node), 2]
  curr <- c(curr, daughters)
  w <- which(x = daughters >= length(x = tree$tip))
  if (length(x = w) > 0) {
    for (i in 1:length(x = w)) {
      curr <- .GetDescendants(tree = tree, node = daughters[w[i]], curr = curr)
    }
  }
  return(curr)
}

#' Function to get all the descendants on a tree left of a given node
#'
#' @param tree Tree object (from ape package)
#' @param node Internal node in the tree
#'
#' @return Returns all descendants left of the given node
#' @noRd
#'
.GetLeftDescendantsSeurat <- function(tree, node) {
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  if (daughters[1] <= (tree$Nnode + 1)) {
    return(daughters[1])
  }
  daughter.use <- .GetDescendants(tree, daughters[1])
  daughter.use <- daughter.use[daughter.use <= (tree$Nnode + 1)]
  return(daughter.use)
}

#' Function to get all the descendants on a tree right of a given node
#'
#' @param tree Tree object (from ape package)
#' @param node Internal node in the tree
#'
#' @return Returns all descendants right of the given node
#' @noRd
.GetRightDescendantsSeurat <- function(tree, node) {
  daughters <- tree$edge[which(x = tree$edge[, 1] == node), 2]
  if (daughters[2] <= (tree$Nnode + 1)) {
    return(daughters[2])
  }
  daughter.use <- .GetDescendants(tree = tree, node = daughters[2])
  daughter.use <- daughter.use[daughter.use <= (tree$Nnode + 1)]
  return(daughter.use)
}

#' Helper function for FindMarkers.Seurat and FoldChange.Seurat
#' Convert idents to cells
#'
#' @importFrom methods is
#' @noRd
#'
.IdentsToCells_Seurat <- function(
    object,
    ident.1,
    ident.2,
    cellnames.use
) {
  #
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  } else if ((length(x = ident.1) == 1 && ident.1[1] == 'clustertree') || is(object = ident.1, class2 = 'phylo')) {
    if (is.null(x = ident.2)) {
      stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
    }
    tree <- if (is(object = ident.1, class2 = 'phylo')) {
      ident.1
    } else {
      SeuratObject::Tool(object = object, slot = 'BuildClusterTree')
    }
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
    }
    ident.1 <- tree$tip.label[.GetLeftDescendantsSeurat(tree = tree, node = ident.2)]
    ident.2 <- tree$tip.label[.GetRightDescendantsSeurat(tree = tree, node = ident.2)]
  }
  if (length(x = as.vector(x = ident.1)) > 1 &&
      any(as.character(x = ident.1) %in% cellnames.use)) {
    bad.cells <- cellnames.use[which(x = !as.character(x = ident.1) %in% cellnames.use)]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- SeuratObject::WhichCells(object = object, idents = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% cellnames.use)) {
    bad.cells <- cellnames.use[which(!as.character(x = ident.2) %in% cellnames.use)]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = cellnames.use, y = ident.1)
    } else {
      ident.2 <- SeuratObject::WhichCells(object = object, idents = ident.2)
    }
  }
  return(list(cells.1 = ident.1, cells.2 = ident.2))
}


#' Differential test with PACS on Seurat object
#' @description
#' Run differential test using the similar style as in Seurat
#' `Seurat::FindMarkers()` function.
#' Note 1: This function is modified from Seurat package v5.1.0
#' Note 2: Seurat must be installed to use this function.
#' @param object A `Seurat` object
#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{BuildClusterTree} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
#' @param group.by Regroup cells into a different identity class prior to
#' performing differential expression (see example)
#' @param subset.ident Subset a particular identity class prior to regrouping.
#' Only relevant if group.by is set (see example)
#' @param assay Assay to use in differential expression testing, for PACS, we
#' use `count` as assay
#' @param reduction Reduction to use in differential expression testing - will
#' test for DE on cell embeddings
#' @param latent.vars Latent variables to be controlled for from the meta.data
#' column of the Seurat object
#' @param ... Additional arguments for identifying differential peaks
#'
#' @concept differential_expression
#' @export
#'
FindMarkersPACS <- function(
    object,
    ident.1 = NULL,
    ident.2 = NULL,
    latent.vars = NULL,
    group.by = NULL,
    subset.ident = NULL,
    assay = "count",
    reduction = NULL,
    ...
) {

  ## make sure object is a Seurat object
  if (!is(object, "Seurat")) {
    stop("Object must be a Seurat object")
  }

  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    SeuratObject::Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  if (length(x = ident.1) == 0) {
    stop("At least 1 ident must be specified in `ident.1`")
  }

  # select which data to use
  if (is.null(x = reduction)) {
    assay <- assay %||% SeuratObject::DefaultAssay(object = object)
    data.use <- object[[assay]]
    cellnames.use <-  colnames(x = data.use)
  } else {
    data.use <- object[[reduction]]
    cellnames.use <- rownames(x = data.use)
  }

  cells <- .IdentsToCells_Seurat(
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
    rlang::abort(
      message = "Cells in one or both identity groups are not present in the data requested"
    )
  }

  # fetch latent.vars
  if (!is.null(x = latent.vars)) {
    latent.vars <- SeuratObject::FetchData(
      object = object,
      vars = latent.vars,
      cells = c(cells$cells.1, cells$cells.2)
    )
  }

  # check normalization method
  norm.command <- paste0("NormalizeData.", assay)
  norm.method <- if (norm.command %in% SeuratObject::Command(object = object) &&
                     is.null(x = reduction)) {
    SeuratObject::Command(
      object = object,
      command = norm.command,
      value = "normalization.method"
    )
  } else if (length(x = intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"),
                                  y = SeuratObject::Command(object = object)))) {
    command <- intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"),
                         y = SeuratObject::Command(object = object))[1]
    SeuratObject::Command(
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

#' Perform differential expression testing using a logistic regression framework
#'
#' Constructs a logistic regression model predicting group membership based on a
#' given feature and compares this to a null model with a likelihood ratio test.
#'
#' @param data.use expression matrix
#' @param cells.1 Vector of cells in group 1
#' @param cells2. Vector of cells in group 2
#' @param latent.vars Latent variables to include in model
#' @param verbose Print messages
#'
#' @importFrom stats as.formula
#' @noRd
.PACSTest <- function(
    data.use,
    cells.1,
    cells.2,
    latent.vars = NULL,
    n_cores = NULL,
    verbose = TRUE
) {

  if("group_test" %in% colnames(latent.vars)) {
    colnames(latent.vars)[colnames(latent.vars) == "group_test"] <- "group_t"
  }

  ## assign group labels
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group_test"] <- "Group1"
  group.info[cells.2, "group_test"] <- "Group2"

  ## get data and meta.data
  data.use <- data.use[, rownames(group.info), drop = FALSE]
  latent.vars <- as.data.frame(latent.vars[rownames(group.info), , drop = FALSE])

  ## calculate capturing rate
  ctypes <- unique(group.info[, "group_test"])
  r_by_ct_out <- PICsnATAC::get_r_by_ct_mat_pq(
    cell_type_set = ctypes,
    r_by_c = data.use,
    cell_type_labels = group.info[, "group_test"],
    n_features_per_cell = dim(data.use)[1]
  )

  ## get formula
  fmla_full <- as.formula(object = paste(
    " ~ group_test +",
    paste(colnames(x = latent.vars), collapse = "+")
  ))

  fmla_null <- as.formula(object = paste(
    " ~ ",
    paste(colnames(x = latent.vars), collapse = "+")
  ))

  ## convert group info into factors
  group.info[, "group_test"] <- factor(x = group.info[, "group_test"])
  latent.vars$group_test <- group.info[, "group_test"]

  if (is.null(n_cores)) {
    cat("automatically determining number of cores", "\n")
    n_cores <- future::nbrOfWorkers()
  }

  pacs_out <- pacs_test_logit(covariate_meta.data = latent.vars,
                              formula_full = fmla_full,
                              formula_null = fmla_null,
                              pic_matrix = data.use,
                              cap_rates = r_by_ct_out$q_vec,
                              n_cores = n_cores)

  p_val <- pacs_out$pacs_p_val

  to.return <- data.frame(p_val, row.names = rownames(data.use))
  return(to.return)
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
    .PACSTest(
      data.use = data.use,
      cells.1 = cells.1,
      cells.2 = cells.2,
      latent.vars = latent.vars,
      verbose = verbose
    )
  return(de.results)
}

