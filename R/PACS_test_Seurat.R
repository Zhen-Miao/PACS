## wrapper for Seurat



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
