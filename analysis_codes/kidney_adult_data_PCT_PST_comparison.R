## test for type1 error simulation

## load libraries
library('SummarizedExperiment')
library('dplyr')
library('Matrix')
library('tidyr')
library('parallel')
library('tictoc')
library('logistf')
library("PICsnATAC")
library("PACS")

## load necessary data
load('data_for_test_for_t1e_power.rdata')  ## kidney data set

## pmats, pmatsbin, unimap2, x.sp_sample2, x.sp_cluster2,ctypes,peak.name

## get labels
cell_type_labels <- x.sp_cluster2
cell_type_set <- c('PT','PT2')


p_by_c <- t(pmatsbin) ## 300755  14458
rm(pmats)
rm(pmatsbin)
n_cell_acc <- Matrix::rowSums(p_by_c)

## store peak names in the matrix
rownames(p_by_c) <- peak.name

quantile(n_cell_acc)
# 0%   25%   50%   75%  100%
# 0    15    45   173 13006

p_by_c <- p_by_c[n_cell_acc >= 15,] ## 229852  14458

# r_by_ct_low = get_r_by_ct_mat_pq(cell_type_set = cell_type_set,
#                                  r_by_c = p_by_c,
#                                  cell_type_labels = cell_type_labels,
#                                  n_features_per_cell = dim(p_by_c)[1],
#                                  p_acc = 0.0005,
#                                  q_acc = 0.0005,
#                                  n_max_iter = 400)
#
#
# saveRDS(r_by_ct_low, 'r_by_ct_low_kidney_2023.rds')
r_by_ct_low = readRDS('r_by_ct_low_kidney_2023.rds')

## keep only relevant cell types
p_by_c = p_by_c[,cell_type_labels %in% cell_type_set] ## 14228 cells
x.sp_sample2 = x.sp_sample2[cell_type_labels %in% cell_type_set] ## 14228 cells
cell_type_labels = cell_type_labels[cell_type_labels %in% cell_type_set] ## 14228 cells

cap_rates = r_by_ct_low$q_vec_new[colnames(p_by_c)]

meta_use = data.frame(id = 1:length(cell_type_labels),
                      CellType = cell_type_labels, sample = x.sp_sample2)

## get dummy variable
xdummy <- meta_use %>% mutate(dummy=1) %>%
  spread(key=sample,value=dummy, fill=0) %>% select(-'90025')
xdummy <- xdummy %>% mutate(dummy=1) %>%
  spread(key=CellType,value=dummy, fill=0) %>% select(-PT)

xdummy <- xdummy %>% select(-id)
xdummy <- as.matrix(xdummy)
xdummy <- cbind(1,xdummy)
colnames(xdummy)[1] <- 'intercept'
x_full <- xdummy


par_initial_one = c(rep(0.05, length = dim(x_full)[2] - 1), c(0.05))


p_by_c = as.matrix(p_by_c)
gc()

full_para = estimate_parameters(r_by_c = p_by_c,
                                design_mat = x_full,
                                par_initial = par_initial_one,
                                cap_rate_vec = cap_rates,
                                mc.cores = 1)

saveRDS(full_para, 'full_para_kidney_adult_PT_PT2.rds')


par_initial_one = c(rep(0.05, length = dim(x_full)[2] - 1), c(0))

fartial_para = estimate_parameters_null(r_by_c = p_by_c,
                                     design_mat = x_full,
                                     par_initial = par_initial_one,
                                     hold_zero = c(4L),
                                     cap_rate_vec = cap_rates,mc.cores = 1)

saveRDS(fartial_para, 'fartial_para_kidney_adult_PT_PT2.rds')


