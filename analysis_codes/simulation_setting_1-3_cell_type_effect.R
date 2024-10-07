## test for type1 error simulation


## load libraries
library('SummarizedExperiment')
library('Matrix')
library('parallel')
library('tictoc')
library('PICsnATAC')
library('PACS')
library(metap)


## load necessary data
load('archR_data_p_q_saved.RData')
lam_list = readRDS('lam_list_low_23.rds') ## list of insertion rate lambda for each cell type

p_by_t_new = do.call(cbind, lam_list)

ct_choice1 <- 'Jurkat'
ct_choice2 <- '293T'

## set n cells
n_cell_total = 4000
n_repeat <- 5


################## choose one of the three settings #####################
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




n_features_per_cell <- dim(p_by_t_new)[1] ## here we use all features for capturing rate estimation
n_features_sample <- 6000 ## we sample some features for testing

########################################
## store the data in the following matrix
#######################################
methods_all <- c('PACS', 'Seurat','archR_n','archR_s',
                 'snapATAC_n','snapATAC_s' ,'Fisher_n', 'Fisher_s')

p_value_mat <- rep(list(), length = length(methods_all))
names(p_value_mat) <- methods_all

for(m in methods_all){
  p_value_mat[[m]] = matrix(ncol = n_repeat, nrow = n_features_sample)
}
conv_mat_our_sub1 <- matrix(ncol = n_repeat, nrow = n_features_sample)




########################################
## get capturing rate
#######################################
q_vec_pos <- q_vec_new[cell_type_labels_n == ct_choice1]
q_vec_neg <- q_vec_new[cell_type_labels_n == ct_choice2]

## filtering -- cells with very small q value should be filtered
q_vec_pos <- q_vec_pos[q_vec_pos > 0.1]
q_vec_neg <- q_vec_neg[q_vec_neg > 0.1]

true_q_total_t1b1 <- sample(q_vec_pos, replace = T, size = n_cell_total)
true_q_total_t1b2 <- sample(q_vec_pos, replace = T, size = n_cell_total)
true_q_total_t2b1 <- sample(q_vec_neg, replace = T, size = n_cell_total)
true_q_total_t2b2 <- sample(q_vec_neg, replace = T, size = n_cell_total)
names(true_q_total_t1b1) <- names(true_q_total_t1b2) <- c(1:n_cell_total)
names(true_q_total_t2b1) <- names(true_q_total_t2b2) <- c(1:n_cell_total)


true_q_total_t1b1[is.na(true_q_total_t1b1)] = 0.5
true_q_total_t1b2[is.na(true_q_total_t1b2)] = 0.5
true_q_total_t2b1[is.na(true_q_total_t2b1)] = 0.5
true_q_total_t2b2[is.na(true_q_total_t2b2)] = 0.5


########################################
## function to get data matrix
#######################################
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

########################################
## get open probability
#######################################
## modify the true open probability
p_by_t_ct1 <- p_by_t_new[,ct_choice1]
p_by_t_ct2 <- p_by_t_new[,ct_choice2]

## remove the most extreme case
p_by_t_ct1[p_by_t_ct1 < 0.14] <- 0.14
p_by_t_ct1[p_by_t_ct1 > 2.5] <- 2.5

p_by_t_ct2[p_by_t_ct2 < 0.14] <- 0.14
p_by_t_ct2[p_by_t_ct2 > 2.5] <- 2.5

mean_p_by_t <- rowMeans(cbind(p_by_t_ct1, p_by_t_ct2))

no_DAR <- which(abs(p_by_t_ct1 - p_by_t_ct2) <= 0.1)
true_DAR <- which(!(abs(p_by_t_ct1 - p_by_t_ct2) <= 0.1))

p_by_t_ct1_mod <- p_by_t_ct1
p_by_t_ct1_mod[no_DAR] <- mean_p_by_t[no_DAR]
p_by_t_ct2_mod <- p_by_t_ct2
p_by_t_ct2_mod[no_DAR] <- mean_p_by_t[no_DAR]

## add batch effect
sample_true_DAR = sample(1:length(true_DAR), replace = T, size = n_features_per_cell/2)
sample_no_DAR = sample(1:length(no_DAR), replace = T, size = n_features_per_cell/2)

p_by_t_ct1_mod_DAR = p_by_t_ct1_mod[true_DAR][sample_true_DAR]
p_by_t_ct1_mod_nDAR = p_by_t_ct1_mod[no_DAR][sample_no_DAR]

p_by_t_ct2_mod_DAR = p_by_t_ct2_mod[true_DAR][sample_true_DAR]
p_by_t_ct2_mod_nDAR = p_by_t_ct2_mod[no_DAR][sample_no_DAR]

## let us try this effect size for now, and try diff ones later on
pt_t1b1 = c(p_by_t_ct1_mod_nDAR, p_by_t_ct1_mod_DAR)
pt_t1b2 = pt_t1b1 * c(rep(1, 23000), rep(0.6,1000),
                                  rep(1.65, 1000),rep(1, 23000),
                              rep(0.6,1000), rep(1.65, 1000) )
pt_t2b1 <- c(p_by_t_ct2_mod_nDAR, p_by_t_ct2_mod_DAR)
pt_t2b2 = pt_t2b1 * c(rep(1, 23000), rep(0.6,1000),
                              rep(1.65, 1000),rep(1, 23000),
                              rep(0.6,1000), rep(1.65, 1000) )

# length(no_DAR) ## 100788
# length(true_DAR) ## 65354

## get the data matrix
## directly sample features here

data_matrix_t1b1_total <- get_data_matrix(true_p_vec = pt_t1b1, true_q_vec = true_q_total_t1b1)
data_matrix_t1b2_total <- get_data_matrix(true_p_vec = pt_t1b2, true_q_vec = true_q_total_t1b2)
data_matrix_t2b1_total <- get_data_matrix(true_p_vec = pt_t2b1, true_q_vec = true_q_total_t2b1)
data_matrix_t2b2_total <- get_data_matrix(true_p_vec = pt_t2b2, true_q_vec = true_q_total_t2b2)
colnames(data_matrix_t1b1_total) <- colnames(data_matrix_t2b1_total) <- c(1:n_cell_total)
colnames(data_matrix_t1b2_total) <- colnames(data_matrix_t2b2_total) <- c(1:n_cell_total)


cells_sampled_mat_t1b1 <- matrix(ncol = n_repeat, nrow = n_cell_sample_t1b1 )
cells_sampled_mat_t1b2 <- matrix(ncol = n_repeat, nrow = n_cell_sample_t1b2 )
cells_sampled_mat_t2b1 <- matrix(ncol = n_repeat, nrow = n_cell_sample_t2b1 )
cells_sampled_mat_t2b2 <- matrix(ncol = n_repeat, nrow = n_cell_sample_t2b2 )
features_sampled <- c(22001:25000, 47001:50000)

for(iii in 1:n_repeat){
  ## random sample some cells
  cells_sampled_mat_t1b1[,iii] <- sample(1:n_cell_total, replace = F,size = n_cell_sample_t1b1)
  cells_sampled_mat_t1b2[,iii] <- sample(1:n_cell_total, replace = F,size = n_cell_sample_t1b2)
  cells_sampled_mat_t2b1[,iii] <- sample(1:n_cell_total, replace = F,size = n_cell_sample_t2b1)
  cells_sampled_mat_t2b2[,iii] <- sample(1:n_cell_total, replace = F,size = n_cell_sample_t2b2)
}


# quantile(abs(p_by_t_ct1_mod - p_by_t_ct2_mod)) ## just to check the code
n_reads_cell_all <- matrix(nrow = 4000, ncol = n_repeat)
## firstly, get the
for(iii in 1:n_repeat){
  ## random sample some cells and peaks
  cells_sampled_t1b1 <- cells_sampled_mat_t1b1[,iii]
  cells_sampled_t1b2 <- cells_sampled_mat_t1b2[,iii]
  cells_sampled_t2b1 <- cells_sampled_mat_t2b1[,iii]
  cells_sampled_t2b2 <- cells_sampled_mat_t2b2[,iii]

  data_matrix_t1b1 <- data_matrix_t1b1_total[,cells_sampled_t1b1]
  data_matrix_t1b2 <- data_matrix_t1b2_total[,cells_sampled_t1b2]
  data_matrix_t2b1 <- data_matrix_t2b1_total[,cells_sampled_t2b1]
  data_matrix_t2b2 <- data_matrix_t2b2_total[,cells_sampled_t2b2]

  n_reads_cell_all[,iii] = c(colSums(data_matrix_t1b1), colSums(data_matrix_t1b2),
                   colSums(data_matrix_t2b1), colSums(data_matrix_t2b2))
}



## then, only keep relevant features
data_matrix_t1b1_total = data_matrix_t1b1_total[features_sampled,]
data_matrix_t1b2_total = data_matrix_t1b2_total[features_sampled,]
data_matrix_t2b1_total = data_matrix_t2b1_total[features_sampled,]
data_matrix_t2b2_total = data_matrix_t2b2_total[features_sampled,]

#----------- input --------------
for(iii in 1:n_repeat){
  ## random sample some cells and peaks
  data_matrix_t1b1 <- data_matrix_t1b1_total[,cells_sampled_mat_t1b1[,iii]]
  data_matrix_t1b2 <- data_matrix_t1b2_total[,cells_sampled_mat_t1b2[,iii]]
  data_matrix_t2b1 <- data_matrix_t2b1_total[,cells_sampled_mat_t2b1[,iii]]
  data_matrix_t2b2 <- data_matrix_t2b2_total[,cells_sampled_mat_t2b2[,iii]]

  true_q_t1b1 <- true_q_total_t1b1[cells_sampled_mat_t1b1[,iii]]
  true_q_t1b2 <- true_q_total_t1b2[cells_sampled_mat_t1b2[,iii]]
  true_q_t2b1 <- true_q_total_t2b1[cells_sampled_mat_t2b1[,iii]]
  true_q_t2b2 <- true_q_total_t2b2[cells_sampled_mat_t2b2[,iii]]

  colnames(data_matrix_t1b1) <- paste("t1b1", c(1:n_cell_sample_t1b1))
  colnames(data_matrix_t1b2) <- paste("t1b2", c(1:n_cell_sample_t1b2))
  colnames(data_matrix_t2b1) <- paste("t2b1", c(1:n_cell_sample_t2b1))
  colnames(data_matrix_t2b2) <- paste("t2b2", c(1:n_cell_sample_t2b2))

  n_reads_cell = n_reads_cell_all[,iii]
  ## our method
  group.info <- c(rep.int(0,times = n_cell_sample_t1b1 + n_cell_sample_t1b2),
                  rep.int(1,times = n_cell_sample_t2b1 + n_cell_sample_t2b2 ))

  batch.info <- c(rep.int(0,times = n_cell_sample_t1b1),
                  rep.int(1, times = n_cell_sample_t1b2),
                  rep.int(0, times = n_cell_sample_t2b1),
                  rep.int(1, times = n_cell_sample_t2b2))

  meta.data <- data.frame(group = group.info, batch = batch.info)

  data_mat = cbind(data_matrix_t1b1, data_matrix_t1b2,
                   data_matrix_t2b1, data_matrix_t2b2)
  data_mat = Matrix(data_mat,sparse = T)
  rownames(data_mat) = paste('f', 1:n_features_sample, sep = '_')

  cap_rates = c(true_q_t1b1, true_q_t1b2,
                true_q_t2b1, true_q_t2b2)

   our_p = pacs_test_sparse(
     covariate_meta.data = meta.data,
     formula_full = ~ factor(group) + factor(batch),
     formula_null = ~ factor(batch),
     pic_matrix = data_mat,
     n_peaks_per_round = NULL,
     T_proportion_cutoff = 0.2,
     cap_rates = cap_rates
   )

  p_value_mat[['PACS']][,iii] = our_p$pacs_p_val


 ## seurat method
 xdummy_null = matrix(rep(1, length(group.info)),ncol = 1)
 seurat_x_dummy_null = cbind(xdummy_null, batch.info, n_reads_cell)
 p_value_mat[['Seurat']][,iii] =
   seurat_method3_subsample(data_mat,seurat_x_dummy_null, group.info)

 ## archR method -- marginal
 p_value_mat[['archR_n']][,iii] = archR_method(cbind(data_matrix_t1b1, data_matrix_t1b2),
                                                      cbind(data_matrix_t2b1, data_matrix_t2b2))

 ## archR method -- stratified
 p_arch_s1 = archR_method(data_matrix_t1b1, data_matrix_t2b1)
 p_arch_s2 = archR_method(data_matrix_t1b2, data_matrix_t2b2)
 p_mat_archR_s = cbind(p_arch_s1, p_arch_s2)
 p_mat_archR_s[p_mat_archR_s >= 1] <- 0.9999
 p_value_mat[['archR_s']][,iii] = apply(X = p_mat_archR_s,1, FUN = function(x) sump(p = x)$p)


 ## snapATAC_marginal
 p_value_mat[['snapATAC_n']][,iii] = snapATAC_method(cbind(data_matrix_t1b1, data_matrix_t1b2),
                                                     cbind(data_matrix_t2b1, data_matrix_t2b2), bcv = 0.4)

 ## snapATAC method -- stratified
 p_sn_s1 = snapATAC_method(data_matrix_t1b1, data_matrix_t2b1, bcv = 0.4)
 p_sn_s2 = snapATAC_method(data_matrix_t1b2, data_matrix_t2b2, bcv = 0.4)
 p_mat_snap_s = cbind(p_sn_s1, p_sn_s2)
 p_mat_snap_s[p_mat_snap_s >= 1] <- 0.9999
 p_value_mat[['snapATAC_s']][,iii] = apply(X = p_mat_snap_s,1, FUN = function(x) sump(p = x)$p)


 ## Fisher
 p_value_mat[['Fisher_n']][,iii] = fisher_method(cbind(data_matrix_t1b1, data_matrix_t1b2),
                                                 cbind(data_matrix_t2b1, data_matrix_t2b2))

 ## Fisher method -- stratified
 p_fi_s1 = fisher_method(data_matrix_t1b1, data_matrix_t2b1)
 p_fi_s2 = fisher_method(data_matrix_t1b2, data_matrix_t2b2)
 p_mat_fi_s = cbind(p_fi_s1, p_fi_s2)
 p_mat_fi_s[p_mat_fi_s == 1] = 0.9999
 p_value_mat[['Fisher_s']][,iii] = apply(X = p_mat_fi_s,1, FUN = function(x) sump(p = x)$p)

}


saveRDS(p_value_mat, 'multifactor_simulation_marginal_p_value_mat_balanced_updated Dec 10 1.65.rds')



## store the data in the following matrix
## store the data in the following matrix
methods_all <- c('PACS', 'Seurat','archR_n','archR_s',
                 'snapATAC_n','snapATAC_s' ,'Fisher_n', 'Fisher_s')


t1power_mat <- matrix(nrow = length(methods_all), ncol = 5)
colnames(t1power_mat) <- c('t1e','t1e_sd','power','power_sd','scenario')
rownames(t1power_mat) <- methods_all
t1power_mat[,'scenario'] <- 1


n_features_sample = 6000

for(jj in methods_all){
  t1e_v <- colMeans(p_value_mat[[jj]][1:(n_features_sample/2),] < 0.05, na.rm = T)
  t1power_mat[jj,'t1e'] <- mean(t1e_v,na.rm = T)
  t1power_mat[jj,'t1e_sd'] <- sd(t1e_v,na.rm = T)

  power_v <- colMeans(p_value_mat[[jj]][(n_features_sample/2+1):n_features_sample,] < 0.05, na.rm = T)
  t1power_mat[jj,'power'] <- mean(power_v,na.rm = T)
  t1power_mat[jj,'power_sd'] <- sd(power_v,na.rm = T)
}

t1power_mat

saveRDS(t1power_mat,'t1power_mat_marginal_Dec10_1.65.rds')



