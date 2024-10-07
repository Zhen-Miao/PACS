
## load libraries
library('SummarizedExperiment')
library('dplyr')
library('Matrix')
library('tidyr')
library('parallel')
library('tictoc')
library('ggplot2')
library('PACS')

meta_use = readRDS('iri_subset_meta_use.rds')
iri = readRDS('iri_subset_pmat.rds')


## get dummy variable

mt = as.character(meta_use$time)
mt[mt == 'sham'] <- 0
mt[mt == '4hr'] <- 1
mt[mt == '12hr'] <- 2


xdummy <- cbind(1L,as.numeric(mt))
xdummy <- as.matrix(xdummy)
colnames(xdummy) <- c('intercept','time')


r_by_ct = get_r_by_ct_mat_pq(cell_type_set = c('ct1'),
                                 r_by_c = iri,
                                 cell_type_labels = rep('ct1', times = ncol(iri)),
                                 n_features_per_cell = dim(iri)[1],
                                 p_acc = 0.0005,
                                 q_acc = 0.0005,
                                 n_max_iter = 400)
saveRDS(r_by_ct,'iri_r_by_ct_PT_subset.rds')

q_vec_new = r_by_ct$q_vec_new[colnames(iri)]

### initialize with random small values -- this step does not matter
par_initial_one = rep(0.02, length = dim(xdummy)[2])

iri = as.matrix(iri)
gc()

full_para = estimate_parameters(r_by_c = iri,
               design_mat = xdummy,
               par_initial = par_initial_one,
               cap_rate_vec = q_vec_new,mc.cores = 1)

saveRDS(full_para, 'full_para_iri.rds')

par_initial_one = c(rep(0.02, length = dim(xdummy)[2]-1), c(0))

partial_para = estimate_parameters_null(r_by_c = iri,
                                        design_mat = xdummy,
                                        par_initial = par_initial_one,
                                        hold_zero = c(2L),
                                        cap_rate_vec = q_vec_new,
                                        mc.cores = 1)

saveRDS(partial_para, 'partial_para_iri.rds')


