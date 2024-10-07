## test for type1 error simulation

## load libraries
library('SummarizedExperiment')
library('dplyr')
library('Matrix')
library('tidyr')
library('parallel')
library('tictoc')
library('PICsnATAC')
library('PACS')

source('code/other_methods_for_differential_updated_April2024.R')



### next step, estimate insertion rate lambda based on this

load('data/archR_data_p_q_saved.RData')
lam_list = readRDS('data/lam_list_low_23.rds')

p_by_t_new = do.call(cbind, lam_list)



ct_choice1 <- 'Jurkat'
ct_choice2 <- '293T'

## set n cells
n_cell_total = 5000
n_repeat <- 5
n_features_per_cell <- 8000 ## here we use all features for capturing rate estimation
n_features_sample <- 5000 ## we sample the same number of t1e and power features
n_features_total <- dim(p_by_t_new)[1]

# q_vec = r_by_ct_low$q_vec_new

# ## balanced cells in group
# n_cell_sample_t1b1 <- 1000
# n_cell_sample_t1b2 <- 1000
# n_cell_sample_t2b1 <- 1000
# n_cell_sample_t2b2 <- 1000

## imbalance cells in group
n_cell_sample_t1b1 <- 1600
n_cell_sample_t1b2 <- 400
n_cell_sample_t2b1 <- 800
n_cell_sample_t2b2 <- 1200
#
# ## very imbalance cells in group
# n_cell_sample_t1b1 <- 1800
# n_cell_sample_t1b2 <- 200
# n_cell_sample_t2b1 <- 400
# n_cell_sample_t2b2 <- 1600


## store the data in the following matrix
methods_all <- c('PACS', 'Seurat', 'Logistic', 'archR_n','archR_s', 'edgeR_pseudo_n', 'edgeR_n',
                 'snapATAC_n','snapATAC_s' ,'Fisher_n', 'Fisher_s')
p_value_mat <- rep(list(), length = length(methods_all))
names(p_value_mat) <- methods_all

for(m in methods_all){
  p_value_mat[[m]] = matrix(ncol = n_repeat, nrow = n_features_sample*2)
}


q_vec_pos <- q_vec_new[cell_type_labels_n == ct_choice1]
q_vec_neg <- q_vec_new[cell_type_labels_n == ct_choice2]

## filtering -- cells with very small q value should be filtered
q_vec_pos <- q_vec_pos[q_vec_pos > 0.1]
q_vec_neg <- q_vec_neg[q_vec_neg > 0.1]



## modify the true open probability
p_by_t_ct1 <- p_by_t_new[,ct_choice1]
p_by_t_ct2 <- p_by_t_new[,ct_choice2]

## remove the most extreme case
p_by_t_ct1[p_by_t_ct1 < 0.1] <- 0.1
p_by_t_ct1[p_by_t_ct1 > 3] <- 3

p_by_t_ct2[p_by_t_ct2 < 0.1] <- 0.1
p_by_t_ct2[p_by_t_ct2 > 3] <- 3

mean_p_by_t <- rowMeans(p_by_t_new[,c(ct_choice1,ct_choice2)])

quantile(p_by_t_ct1 - p_by_t_ct2)
quantile(abs(p_by_t_ct1 - p_by_t_ct2))

no_DAR <- which(abs(p_by_t_ct1 - p_by_t_ct2) <= 0.1)
true_DAR <- which(!(abs(p_by_t_ct1 - p_by_t_ct2) <= 0.1))

p_by_t_ct1_mod <- p_by_t_ct1
p_by_t_ct1_mod[no_DAR] <- mean_p_by_t[no_DAR]
p_by_t_ct2_mod <- p_by_t_ct2
p_by_t_ct2_mod[no_DAR] <- mean_p_by_t[no_DAR]

f_sample = sample(1:length(true_DAR), size = n_features_per_cell, replace = F)
f_sample2 = sample(1:length(no_DAR), size = n_features_per_cell, replace = F)


p_by_t_ct1_mod_DAR = p_by_t_ct1_mod[true_DAR][f_sample]
p_by_t_ct1_mod_nDAR = p_by_t_ct1_mod[no_DAR][f_sample2]

p_by_t_ct2_mod_DAR = p_by_t_ct2_mod[true_DAR][f_sample]
p_by_t_ct2_mod_nDAR = p_by_t_ct2_mod[no_DAR][f_sample2]

p_by_t_ct1_mod_sorted <- c(p_by_t_ct1_mod_nDAR, p_by_t_ct1_mod_DAR)
p_by_t_ct2_mod_sorted <- c(p_by_t_ct2_mod_nDAR, p_by_t_ct2_mod_DAR)

# length(no_DAR) ## 100788  --> 18026
# length(true_DAR) ## 65354  --> 31974
## store the data matrix

get_data_matrix <- function(true_p_vec,
                            true_q_vec
){
  n_features_per_cell <-  length(true_p_vec)
  n_cell_total <-  length(true_q_vec)
  data_matrix <- matrix(nrow = n_features_per_cell, ncol = n_cell_total)

  for(ii in 1:n_features_per_cell){
    dist = .get_theoretical_c1(insertion_rate = true_p_vec[ii])
    dist2 = c(dist[1],dist[2], sum(dist[3:6]))
    s = rmultinom(n = n_cell_total, size = 1, prob = dist2)
    s[2,] = s[2,] * 2
    s[3,] = s[3,] * 3
    data_matrix[ii,] = colSums(s) - 1
  }

  capturing_matrix <- matrix(nrow = n_features_per_cell, ncol = n_cell_total)
  for(ss in 1:n_cell_total){
    capturing_matrix[,ss] = rbinom(n = n_features_per_cell, size = 1, prob = true_q_vec[ss])
  }

  data_matrix = data_matrix * capturing_matrix
  return(data_matrix)
}


true_q_pos_total <- sample(q_vec_pos, replace = T, size = n_cell_total)
true_q_neg_total <- sample(q_vec_neg, replace = T, size = n_cell_total)
names(true_q_pos_total) <- names(true_q_neg_total) <- c(1:n_cell_total)

true_q_pos_total[is.na(true_q_pos_total)] = 0.5
true_q_neg_total[is.na(true_q_neg_total)] = 0.5

## get the data matrix
data_matrix_pos_total <- get_data_matrix(true_p_vec = p_by_t_ct1_mod_sorted,
                                         true_q_vec = true_q_pos_total)
data_matrix_neg_total <- get_data_matrix(true_p_vec = p_by_t_ct2_mod_sorted,
                                         true_q_vec = true_q_neg_total)
colnames(data_matrix_pos_total) <- colnames(data_matrix_neg_total) <- c(1:n_cell_total)


cells_sampled1_mat <- matrix(ncol = n_repeat, nrow = n_cell_sample_t1b1 )
cells_sampled2_mat <- matrix(ncol = n_repeat, nrow = n_cell_sample_t1b2  )
cells_sampled3_mat <- matrix(ncol = n_repeat, nrow = n_cell_sample_t2b1 )
cells_sampled4_mat <- matrix(ncol = n_repeat, nrow = n_cell_sample_t2b2 )
features_sampled_mat <- matrix(ncol = n_repeat, nrow = n_features_sample *2 )

for(iii in 1:5){
  ## random sample some cells
  cells_sampled1_mat[,iii] <- sample(1:n_cell_total, replace = F,size = n_cell_sample_t1b1)
  cells_sampled2_mat[,iii] <- sample(1:n_cell_total, replace = F,size = n_cell_sample_t1b2)
  cells_sampled3_mat[,iii] <- sample(1:n_cell_total, replace = F,size = n_cell_sample_t2b1)
  cells_sampled4_mat[,iii] <- sample(1:n_cell_total, replace = F,size = n_cell_sample_t2b2)
  f1 = sample(1:n_features_per_cell, replace = F, size = n_features_sample)
  f2 = sample(1:n_features_per_cell, replace = F, size = n_features_sample) +
    n_features_per_cell
  features_sampled_mat[,iii] <- c(f1,f2)
}



# quantile(abs(p_by_t_ct1_mod - p_by_t_ct2_mod)) ## just to check the code



#----------- input --------------
for(iii in 1:5){
  ## random sample some cells and peaks
  cells_sampled1 <- cells_sampled1_mat[,iii]
  cells_sampled2 <- cells_sampled2_mat[,iii]
  cells_sampled3 <- cells_sampled3_mat[,iii]
  cells_sampled4 <- cells_sampled4_mat[,iii]

  ### the two groups are reversed, so that we need to test the
  ### interaction effect
  data_matrix_1 <- data_matrix_pos_total[,cells_sampled1]
  data_matrix_2 <- data_matrix_neg_total[,cells_sampled2]
  data_matrix_3 <- data_matrix_neg_total[,cells_sampled3]

  ## half of features have interaction effect of a-b-b-a,
  ## the other half have a-b-b-b
  data_matrix_4 <- data_matrix_pos_total[,cells_sampled4]
  data_matrix_4_2 <- data_matrix_neg_total[,cells_sampled4]
  data_matrix_4[c(1:n_features_per_cell/2,
                  (n_features_per_cell+1):
                    (n_features_per_cell*1.5)),] <-
    data_matrix_4_2[c(1:n_features_per_cell/2,
                    (n_features_per_cell+1):(n_features_per_cell*1.5)),]

  true_q_1 <- true_q_pos_total[cells_sampled1]
  true_q_2 <- true_q_neg_total[cells_sampled2]
  true_q_3 <- true_q_neg_total[cells_sampled3]
  true_q_4 <- true_q_pos_total[cells_sampled4]

  # print(head(cells_sampled1))

  colnames(data_matrix_1) <- c(1:n_cell_sample_t1b1)
  colnames(data_matrix_2) <- c(1:n_cell_sample_t1b2)
  colnames(data_matrix_3)  <- c(1:n_cell_sample_t2b1)
  colnames(data_matrix_4) <- c(1:n_cell_sample_t2b2)

  ## get some quantities
  n_1 = dim(data_matrix_1)[2]
  n_2 = dim(data_matrix_2)[2]
  n_3 = dim(data_matrix_3)[2]
  n_4 = dim(data_matrix_4)[2]
  n_reads_cell = c(colSums(data_matrix_1), colSums(data_matrix_2),
                   colSums(data_matrix_3), colSums(data_matrix_4))

  ## sample features
  data_matrix_1 = data_matrix_1[features_sampled_mat[,iii],]
  data_matrix_2 = data_matrix_2[features_sampled_mat[,iii],]
  data_matrix_3 = data_matrix_3[features_sampled_mat[,iii],]
  data_matrix_4 = data_matrix_4[features_sampled_mat[,iii],]

  # ## binarize
  # data_matrix_pos_bin <- ifelse(data_matrix_pos != 0, 1, 0)
  # data_matrix_neg_bin <- ifelse(data_matrix_neg != 0, 1, 0)

  ## our method
  group.info <- c(rep.int(0,times = n_1 + n_2),
                  rep.int(1,times = n_3 + n_4 ))

  batch.info <- c(rep.int(0,times = n_1),
                  rep.int(1, times = n_2),
                  rep.int(0, times = n_3),
                  rep.int(1, times = n_4))

  meta.data <- data.frame(group = group.info, batch = batch.info)

  data_mat = cbind(data_matrix_1, data_matrix_2,data_matrix_3, data_matrix_4)
  data_mat = Matrix(data_mat,sparse = T)
  rownames(data_mat) = paste('f', 1:(n_features_sample*2), sep = '_')

  cap_rates = c(true_q_1, true_q_2, true_q_3, true_q_4)

  our_p = pacs_test_sparse(
    covariate_meta.data = meta.data,
    formula_full = ~ group * batch ,
    formula_null = ~ group + batch,
    pic_matrix = data_mat,
    n_peaks_per_round = NULL,
    T_proportion_cutoff = 0.2,
    cap_rates = cap_rates
  )
 #
 p_value_mat[['PACS']][,iii] = our_p$pacs_p_val
 #
 ## seurat
 xdummy_null = matrix(rep(1, length(group.info)),ncol = 1)
 seurat_x_dummy_null = cbind(xdummy_null, batch.info,n_reads_cell)
   ss = seurat_method3_subsample(data_mat, seurat_x_dummy_null, group.info)


 p_value_mat[['Seurat']][,iii] = ss

 ## logistic regression
 seurat_x_dummy_null_2 = cbind(batch.info, n_reads_cell)
 p_value_mat[['Logistic']][,iii] =
   standard_logit_subsample(data_mat,seurat_x_dummy_null_2, group.info)

 ## archR method -- naive
 p_value_mat[['archR_n']][,iii] = archR_method(cbind(data_matrix_1, data_matrix_2),
                                                cbind(data_matrix_3, data_matrix_4))

 ## archR method -- stratified
 p_arch_s1 = archR_method(data_matrix_1, data_matrix_3)
 p_arch_s2 = archR_method(data_matrix_2, data_matrix_4)
 p_mat_archR_s = cbind(p_arch_s1, p_arch_s2)
 p_mat_archR_s[p_mat_archR_s >= 1] <- 0.9999

 p_value_mat[['archR_s']][,iii] = apply(X = p_mat_archR_s,1, FUN = function(x) sumz(p = x)$p)

 ## snapATAC_marginal
 p_value_mat[['snapATAC_n']][,iii] = snapATAC_method(cbind(data_matrix_1, data_matrix_2),
                                                     cbind(data_matrix_3, data_matrix_4), bcv = 0.4)
 ## snapATAC method -- stratified
 p_sn_s1 = snapATAC_method(data_matrix_1, data_matrix_3, bcv = 0.4)
 p_sn_s2 = snapATAC_method(data_matrix_2, data_matrix_4, bcv = 0.4)
 p_mat_snap_s = cbind(p_sn_s1, p_sn_s2)
 p_mat_snap_s[p_mat_snap_s >= 1] <- 0.9999
 p_value_mat[['snapATAC_s']][,iii] = apply(X = p_mat_snap_s,1, FUN = function(x) sumz(p = x)$p)

 #######################################
 ## edgeR pseudobulk  -- considering the complex study design
 #######################################

 mat_combined = cbind(data_matrix_1, data_matrix_2,
                      data_matrix_3, data_matrix_4)
 ## edgeR-multi method
 group.info <- c(rep.int(0,times = n_cell_sample_t1b1 + n_cell_sample_t1b2),
                 rep.int(1,times = n_cell_sample_t2b1 + n_cell_sample_t2b2 ))
 batch.info <- c(rep.int(0,times = n_cell_sample_t1b1),
                 rep.int(1, times = n_cell_sample_t1b2),
                 rep.int(0, times = n_cell_sample_t2b1),
                 rep.int(1, times = n_cell_sample_t2b2))
 meta.data <- data.frame(group = group.info, batch = batch.info)


 p_value_mat[['edgeR_n']][,iii] = edgeR_multi_group(matrix_combined = mat_combined,
                                                    design_mat = meta.data,
                                                    coef_test = 'group')

 pseudo_t1_b1 <- create_pseudo_bulk(data_matrix_1, 10)
 pseudo_t1_b2 <- create_pseudo_bulk(data_matrix_2, 10)
 pseudo_t2_b1 <- create_pseudo_bulk(data_matrix_3, 10)
 pseudo_t2_b2 <- create_pseudo_bulk(data_matrix_4, 10)
 pseudo_combined = cbind(pseudo_t1_b1, pseudo_t1_b2,
                         pseudo_t2_b1, pseudo_t2_b2)

 group.info_pseudo <- c(rep.int(0,times = ncol(pseudo_t1_b1) + ncol(pseudo_t1_b2)),
                        rep.int(1,times = ncol(pseudo_t2_b1) + ncol(pseudo_t2_b2)))
 batch.info_pseudo <- c(rep.int(0,times = ncol(pseudo_t1_b1)),
                        rep.int(1, times = ncol(pseudo_t1_b2)),
                        rep.int(0, times = ncol(pseudo_t2_b1)),
                        rep.int(1, times = ncol(pseudo_t2_b2)))
 meta.data_pseudo <- data.frame(group = group.info_pseudo, batch = batch.info_pseudo)

 p_value_mat[['edgeR_pseudo_n']][,iii] = edgeR_multi_group(
   matrix_combined = pseudo_combined,
   design_mat = meta.data_pseudo,
   coef_test = 'group'
 )

 ## Fisher
 p_value_mat[['Fisher_n']][,iii] = fisher_method(cbind(data_matrix_1, data_matrix_2),
                                                 cbind(data_matrix_3, data_matrix_4))

 ## Fisher method -- stratified
 p_fi_s1 = fisher_method(data_matrix_1, data_matrix_3)
 p_fi_s2 = fisher_method(data_matrix_2, data_matrix_4)
 p_mat_fisher_s = cbind(p_fi_s1, p_fi_s2)
 p_mat_fisher_s[p_mat_fisher_s >= 1] <- 0.9999
 p_value_mat[['Fisher_s']][,iii] = apply(X = p_mat_fisher_s,1, FUN = function(x) sumz(p = x)$p)

}


saveRDS(p_value_mat, 'results/interaction_effect_simulation_unbalanced Apr30.rds')


