
## load libraries
library('edgeR')
library('presto')
library('SummarizedExperiment')
library('dplyr')
library('aricode')
library('Matrix')
library('PICsnATAC')
library('PACS')
library('Rfast')
library(viridisLite);
library(ggplot2);
library("RColorBrewer")
library(dplyr)
library(reshape2)
library(aricode)
library('GenomicRanges')
library(reshape2)


x.sp2 = readRDS('all_comb snapATAC after batch correction with peak info new clu.RDS')

## get rid of P0 data -- 90028 and 90029
x.sp2@cluster[!(x.sp2@sample %in% c('90028', '90029'))] <- 'na'
ctypes <- c('NP','PT',"PC","DCT","immune",'Endo','PT2','IC','Podo','LOH','stroma2')
peak.name <- x.sp2@peak$name


pmats = readRDS('pmat of KS data before binarization.rds')
pmats= pmats[x.sp2@sample %in% c('90028', '90029') & x.sp2@cluster %in% ctypes,]
unimap2 = as.numeric(x.sp2@metaData$UQ[x.sp2@sample %in% c('90028', '90029') &
                                         x.sp2@cluster %in% ctypes])

x.sp_sample2 = x.sp2@sample[x.sp2@sample %in% c('90028', '90029')& x.sp2@cluster %in% ctypes]
x.sp_cluster2 = as.character(x.sp2@cluster[x.sp2@sample %in% c('90028', '90029') &
                                             x.sp2@cluster %in% ctypes])

pmatsbin <- pmats
pmatsbin@x <- rep(1, times = length(pmatsbin@x))


n_cell_acc <- Matrix::colSums(pmatsbin)
pmatsbin <- pmatsbin[,n_cell_acc >= 20]

kidney_mat <- t(pmatsbin)

cell_type_labels <- x.sp_cluster2
cell_type_set <- unique(x.sp_cluster2)
pks = x.sp2@peak$name[n_cell_acc >= 20]
rownames(kidney_mat) <- pks


kidney_mat = readRDS('kidney_mat_P0_2023.rds')
cell_type_labels = readRDS('kidney_cell_type_labels_P0_2023.rds')
cell_type_set = unique(cell_type_labels)

##################################
## stratified sampling
##################################
for(i in 1:10){
  cell_label_all = data.frame(cells = colnames(kidney_mat), labels = cell_type_labels)

  stratified <- cell_label_all %>%
    group_by(labels) %>%
    sample_frac(size = .1)

  testing_set = stratified$cells
  training_set = setdiff(colnames(kidney_mat), testing_set)

  testing_sel <- colnames(kidney_mat) %in% testing_set
  training_sel <- colnames(kidney_mat) %in% training_set



  r_by_ct_low = get_r_by_ct_mat_pq(cell_type_set = cell_type_set,
                                   r_by_c = kidney_mat[,training_sel],
                                   cell_type_labels = cell_type_labels[training_sel],
                                   n_features_per_cell = dim(kidney_mat)[1],
                                   p_acc = 0.0005,
                                   q_acc = 0.0005,
                                   n_max_iter = 400)

  ## identify variable peaks by standard deviation
  p_by_t_low = as.matrix(r_by_ct_low$p_by_t_new)
  p_sd = sqrt(Rfast::rowVars(p_by_t_low))
  quantile(p_sd)


  pk_sel = p_sd > 0.02  ## 50%


  est_low_low = estimate_label_no_cap_rate(r_by_t = r_by_ct_low$p_by_t_new[pk_sel,],
                                           in_r_by_c = kidney_mat[pk_sel,testing_sel] )
  # saveRDS(est_low_low,'est_kidney_mat_9CV_round1.rds')

  ## or only look at singlets
  esti_m5 <- est_low_low

  ### check the overall distribution
  edf = as.data.frame(esti_m5)
  edf$labelss = cell_type_labels[testing_sel]
  ctypes <- colnames(esti_m5)


  ## estimate cell type labels based on the updated version
  elems1 <- rep(1, length = dim(esti_m5)[1])
  elems2 <- rep(2, length = dim(esti_m5)[1])
  esti_ctype_fir = rownth(esti_m5,elems1, num.of.nths = 1,descending = TRUE,index.return = TRUE, parallel = FALSE)
  esti_ctype_sec = rownth(esti_m5,elems2, num.of.nths = 1,descending = TRUE,index.return = TRUE, parallel = FALSE)

  true_ctype_labels = cell_type_labels[testing_sel]

  fst_esti_ctype_labels = ctypes[esti_ctype_fir]
  snd_esti_ctype_labels = ctypes[esti_ctype_sec]

  # ARI(fst_esti_ctype_labels, true_ctype_labels)
  # table(fst_esti_ctype_labels, true_ctype_labels)

  ari = ARI(fst_esti_ctype_labels, true_ctype_labels)
  print('is ARI')
  print(ari)
  print('is ARI')
  print(table(fst_esti_ctype_labels, true_ctype_labels))
  # table(fst_esti_ctype_labels[true_ctype_labels != 'DBL'], true_ctype_labels[true_ctype_labels != 'DBL'])
  print('round')
  print(i)

  print(clustComp(fst_esti_ctype_labels, true_ctype_labels))

}


tb = table(fst_esti_ctype_labels, true_ctype_labels)
write.csv(tb, 'confusion_kidney_P0_PACS.csv')

## check some examples
table(true_ctype_labels, meta$area)



##################################
## broad definition cell lineage
##################################

cell_type_labels <- meta$CellType
cell_type_labels[cell_type_labels == 'ulEN'] <- 'Ex_Neurons'
cell_type_labels[cell_type_labels == 'dlEN'] <- 'Ex_Neurons'
cell_type_labels[cell_type_labels == 'earlyEN'] <- 'Ex_Neurons'
cell_type_labels[cell_type_labels == 'IPC'] <- 'Ex_Neurons'
cell_type_labels[cell_type_labels == 'RG'] <- 'Ex_Neurons'
cell_type_labels[cell_type_labels == 'Insular_Neurons'] <- 'Ex_Neurons'

cell_type_labels[cell_type_labels == 'IN_CGE'] <- 'In_Neurons'
cell_type_labels[cell_type_labels == 'IN_MGE'] <- 'In_Neurons'
cell_type_labels[cell_type_labels == 'MGE_Progenitors'] <- 'In_Neurons'

cell_type_labels[cell_type_labels == 'AstroOligo'] <- 'Astro_Oligo'
cell_type_labels[cell_type_labels == 'Microglia'] <- 'Microglia_'
cell_type_labels[cell_type_labels == 'EndoMural'] <- 'Endo_Mural'

cell_type_set = unique(cell_type_labels)

for(i in 1:10){
  cell_label_all = data.frame(cells = colnames(kidney_mat), labels = cell_type_labels)

  stratified <- cell_label_all %>%
    group_by(labels) %>%
    sample_frac(size = .1)

  testing_set = stratified$cells
  training_set = setdiff(colnames(kidney_mat), testing_set)

  testing_sel <- colnames(kidney_mat) %in% testing_set
  training_sel <- colnames(kidney_mat) %in% training_set



  r_by_ct_low = get_r_by_ct_mat_pq(cell_type_set = cell_type_set,
                                   r_by_c = kidney_mat[,training_sel],
                                   cell_type_labels = cell_type_labels[training_sel],
                                   n_features_per_cell = dim(kidney_mat)[1],
                                   p_acc = 0.0005,
                                   q_acc = 0.0005,
                                   n_max_iter = 400)

  ## identify variable peaks by standard deviation
  p_by_t_low = as.matrix(r_by_ct_low$p_by_t_new)
  p_sd = sqrt(Rfast::rowVars(p_by_t_low))
  quantile(p_sd)


  pk_sel = p_sd > 0.02




  est_low_low = estimate_label_no_cap_rate(r_by_t = r_by_ct_low$p_by_t_new[pk_sel,],
                                           in_r_by_c = kidney_mat[pk_sel,testing_sel] )
  # saveRDS(est_low_low,'est_kidney_mat_9CV_round1.rds')

  ## or only look at singlets
  esti_m5 <- est_low_low

  ### check the overall distribution
  edf = as.data.frame(esti_m5)
  edf$labelss = cell_type_labels[testing_sel]
  ctypes <- colnames(esti_m5)


  ## estimate cell type labels based on the updated version
  elems1 <- rep(1, length = dim(esti_m5)[1])
  elems2 <- rep(2, length = dim(esti_m5)[1])
  esti_ctype_fir = rownth(esti_m5,elems1, num.of.nths = 1,descending = TRUE,index.return = TRUE, parallel = FALSE)
  esti_ctype_sec = rownth(esti_m5,elems2, num.of.nths = 1,descending = TRUE,index.return = TRUE, parallel = FALSE)

  true_ctype_labels = cell_type_labels[testing_sel]

  fst_esti_ctype_labels = ctypes[esti_ctype_fir]
  snd_esti_ctype_labels = ctypes[esti_ctype_sec]

  # ARI(fst_esti_ctype_labels, true_ctype_labels)
  # table(fst_esti_ctype_labels, true_ctype_labels)

  ari = ARI(fst_esti_ctype_labels, true_ctype_labels)
  print('is ARI')
  print(ari)
  print('is ARI')
  print(table(fst_esti_ctype_labels, true_ctype_labels))
  # table(fst_esti_ctype_labels[true_ctype_labels != 'DBL'], true_ctype_labels[true_ctype_labels != 'DBL'])
  print('round')
  print(i)

  print(clustComp(fst_esti_ctype_labels, true_ctype_labels))

}

##################################
## stratified sampling -- low to high
##################################


  r_by_ct_low = get_r_by_ct_mat_pq(cell_type_set = cell_type_set,
                                   r_by_c = kidney_mat,
                                   cell_type_labels = cell_type_labels,
                                   n_features_per_cell = 234171,
                                   p_acc = 0.0005,
                                   q_acc = 0.0005,
                                   n_max_iter = 400)

  ## identify variable peaks by standard deviation
  p_by_t = as.matrix(r_by_ct_low$p_by_t_new)
  p_sd = sqrt(Rfast::rowVars(p_by_t))
  quantile(p_sd)

  pk_sel = p_sd > 0.04

  est_low_h = estimate_label_no_cap_rate(r_by_t = r_by_ct_low_h$p_by_t_new[pk_sel,],
                                       in_r_by_c = kidney_mat_h[pk_sel,] )
  # saveRDS(est_low_low,'est_low_low_2023_9fcv.rds')

  ## or only look at singlets
  esti_m5 <- est_low_h

  ### check the overall distribution
  edf = as.data.frame(esti_m5)
  edf$labelss = cell_type_labels_h
  ctypes <- colnames(esti_m5)


  ## estimate cell type labels based on the updated version
  elems1 <- rep(1, length = dim(esti_m5)[1])
  elems2 <- rep(2, length = dim(esti_m5)[1])
  esti_ctype_fir = rownth(esti_m5,elems1, num.of.nths = 1,descending = TRUE,index.return = TRUE, parallel = FALSE)
  esti_ctype_sec = rownth(esti_m5,elems2, num.of.nths = 1,descending = TRUE,index.return = TRUE, parallel = FALSE)

  true_ctype_labels_h = cell_type_labels_h

  fst_esti_ctype_labels = ctypes[esti_ctype_fir]
  snd_esti_ctype_labels = ctypes[esti_ctype_sec]

  # ARI(fst_esti_ctype_labels, true_ctype_labels)
  # table(fst_esti_ctype_labels, true_ctype_labels)

  ari = ARI(fst_esti_ctype_labels[true_ctype_labels_h != 'DBL'], true_ctype_labels_h[true_ctype_labels_h != 'DBL'])
  print(ari)
  # table(fst_esti_ctype_labels[true_ctype_labels != 'DBL'], true_ctype_labels[true_ctype_labels != 'DBL'])

##################################
## stratified sampling end
##################################


## filtering based on p_ct_filt
rowsd <- rowSds(p_ct_filt)
elems1 <- rep(1, length = dim(p_ct_filt)[1])
high_prop = rownth(p_ct_filt,elems1, num.of.nths = 1,
                        descending = TRUE,index.return = F, parallel = FALSE)
low_prop = rownth(p_ct_filt,elems1, num.of.nths = 1,
                   descending = F,index.return = F, parallel = FALSE)

diff_prop <- high_prop - low_prop
# hist(diff_prop)
kidney_mat = kidney_mat[diff_prop >= 0.05 & rowsd >= 0.015 ,] ## 188055 features left
kidney_mat_h = kidney_mat_h[diff_prop >= 0.05 & rowsd >= 0.015 ,] ## 188055 features left



p_ct_ave <- get_kidney_matt_mat_ave(kidney_mat = kidney_mat,group_set = cell_type_set,group_lables = cell_type_labels,
                                capturing_rate =  q_vec,max_p = 0.999, min_p = 0.001,adding_doublet = T)

save.image(file = "ArchR_whole.RData")

est_m_doub <- estimate_type(p_by_t = p_ct_half, in_kidney_mat = kidney_mat,
                            capturing_rate = q_vec, alpha = 1 )


saveRDS(est_m_doub, 'est_m_ArchR_doub_half.rds')






