

## load libraries
library('presto')
library('SummarizedExperiment')
library('dplyr')
library('aricode')
library('PICsnATAC')
library('PACS')

## the data are from GSE162690
low = readRDS('GSE162690_CellLine_LowLoading-PeakMatrix-SE.rds')
high = readRDS('GSE162690_CellLine_HighLoading-PeakMatrix-SE.rds')


## get the number of uniquely mapped reads
unq_reads <- low@colData@listData$nMonoFrags
all_reads <- low@colData@listData$nFrags

## get the number of uniquely mapped reads
unq_reads_h <- high@colData@listData$nMonoFrags
all_reads_h <- high@colData@listData$nFrags




######################
### load labels and matrix
#######################
## this is the gold annotation
table(low@colData@listData$DemuxletClassify2)  ## AMB 5, DBL 336
cell_type_labels <- low@colData@listData$DemuxletClassify2
cell_type_set <- unique(cell_type_labels)
cell_type_set <- setdiff(cell_type_set,c('AMB','DBL'))

p_by_c <- low@assays$data@listData$PeakMatrix ## 336098  10832
p_by_c@x = rep(1, length = length(p_by_c@x))
n_cell_acc <- Matrix::rowSums(p_by_c)

quantile(n_cell_acc)
# 0%  25%  50%  75% 100%
# 0   37   66  138 3041

##########################
## the same for high reads
## this is the gold annotation
table(high@colData@listData$DemuxletClassify2)  ## AMB 16, DBL 1056
cell_type_labels_h <- high@colData@listData$DemuxletClassify2
cell_type_set_h <- unique(cell_type_labels_h)
cell_type_set_h <- setdiff(cell_type_set_h,c('AMB','DBL'))

p_by_c_h <- high@assays$data@listData$PeakMatrix ## 336098  14966
p_by_c_h@x = rep(1, length = length(p_by_c_h@x))
n_cell_acc_h <- Matrix::rowSums(p_by_c_h)

quantile(n_cell_acc_h)
# 0%  25%  50%  75% 100%
# 1   49   82  164 3860

q1 = quantile(n_cell_acc,.25)
q2 = quantile(n_cell_acc_h,.25)

## get peak names
pks = as.data.frame(low@rowRanges)
pks2 = paste(pks$seqnames,':',pks$start,'-',pks$end,sep = '')
pks2 = pks2[n_cell_acc >= q1 & n_cell_acc_h >= q2]


### filtering peaks
p_by_c <- p_by_c[n_cell_acc >= q1 & n_cell_acc_h >= q2,] ## 234171 peaks  10832 cells
p_by_c_h <- p_by_c_h[n_cell_acc >= q1 & n_cell_acc_h >= q2,] ## 234171 peaks  16047 cells

p_by_c <- p_by_c[,cell_type_labels != 'AMB']    ## 11168
p_by_c_h <- p_by_c_h[,cell_type_labels_h != 'AMB']    ## 16031

rownames(p_by_c) <- pks2
rownames(p_by_c_h) <- pks2

cell_type_labels = cell_type_labels[cell_type_labels != 'AMB']
cell_type_labels_h = cell_type_labels_h[cell_type_labels_h != 'AMB']

##################################
## stratified sampling
##################################
for(i in 1:7){
  cell_label_all = data.frame(cells = colnames(p_by_c), labels = cell_type_labels)

  stratified <- cell_label_all %>%
    group_by(labels) %>%
    sample_frac(size = .1)

  testing_set = stratified$cells
  training_set = setdiff(colnames(p_by_c), testing_set)

  testing_sel <- colnames(p_by_c) %in% testing_set
  training_sel <- colnames(p_by_c) %in% training_set



  r_by_ct_low = get_r_by_ct_mat_pq(cell_type_set = cell_type_set,
                                   r_by_c = p_by_c[,training_sel],
                                   cell_type_labels = cell_type_labels[training_sel],
                                   n_features_per_cell = 234171,
                                   p_acc = 0.0005,
                                   q_acc = 0.0005,
                                   n_max_iter = 400)

  ## identify variable peaks by standard deviation
  p_by_t_low = as.matrix(r_by_ct_low$p_by_t_new)
  p_sd = sqrt(Rfast::rowVars(p_by_t_low))
  quantile(p_sd)

  pk_sel = p_sd > 0.04

  est_low_low = estimate_label_no_cap_rate(r_by_t = r_by_ct_low$p_by_t_new[pk_sel,],
                                           in_r_by_c = p_by_c[pk_sel,testing_sel] )
  # saveRDS(est_low_low,'est_low_low_2023_9fcv.rds')

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

  ari = ARI(fst_esti_ctype_labels[true_ctype_labels != 'DBL'], true_ctype_labels[true_ctype_labels != 'DBL'])
  print(ari)
  table(fst_esti_ctype_labels[true_ctype_labels != 'DBL'], true_ctype_labels[true_ctype_labels != 'DBL'])
  print('round')
  print(i)
}


tb = table(fst_esti_ctype_labels[true_ctype_labels != 'DBL'], true_ctype_labels[true_ctype_labels != 'DBL'])

write.csv(tb, 'confusion_cell_line_low_PACS.csv')

##################################
## stratified sampling end
##################################

##################################
## stratified sampling -- high to high
##################################
for(i in 1:7){
  cell_label_all_h = data.frame(cells = colnames(p_by_c_h), labels = cell_type_labels_h)

  stratified_h <- cell_label_all_h %>%
    group_by(labels) %>%
    sample_frac(size = .1)

  testing_set_h = stratified_h$cells
  training_set_h = setdiff(colnames(p_by_c_h), testing_set_h)

  testing_sel_h <- colnames(p_by_c_h) %in% testing_set_h
  training_sel_h <- colnames(p_by_c_h) %in% training_set_h



  r_by_ct_low_h = get_r_by_ct_mat_pq(cell_type_set = cell_type_set,
                                   r_by_c = p_by_c_h[,training_sel_h],
                                   cell_type_labels = cell_type_labels_h[training_sel_h],
                                   n_features_per_cell = 234171,
                                   p_acc = 0.0005,
                                   q_acc = 0.0005,
                                   n_max_iter = 400)

  ## identify variable peaks by standard deviation
  p_by_t_h = as.matrix(r_by_ct_low_h$p_by_t_new)
  p_sd_h = sqrt(Rfast::rowVars(p_by_t_h))
  quantile(p_sd_h)

  pk_sel_h = p_sd_h > 0.04

  est_h_h = estimate_label_no_cap_rate(r_by_t = r_by_ct_low_h$p_by_t_new[pk_sel_h,],
                                           in_r_by_c = p_by_c_h[pk_sel_h,testing_sel_h] )
  # saveRDS(est_low_low,'est_low_low_2023_9fcv.rds')

  ## or only look at singlets
  esti_m5 <- est_h_h

  ### check the overall distribution
  edf = as.data.frame(esti_m5)
  edf$labelss = cell_type_labels_h[testing_sel_h]
  ctypes <- colnames(esti_m5)


  ## estimate cell type labels based on the updated version
  elems1 <- rep(1, length = dim(esti_m5)[1])
  elems2 <- rep(2, length = dim(esti_m5)[1])
  esti_ctype_fir = rownth(esti_m5,elems1, num.of.nths = 1,descending = TRUE,index.return = TRUE, parallel = FALSE)
  esti_ctype_sec = rownth(esti_m5,elems2, num.of.nths = 1,descending = TRUE,index.return = TRUE, parallel = FALSE)

  true_ctype_labels_h = cell_type_labels_h[testing_sel_h]

  fst_esti_ctype_labels = ctypes[esti_ctype_fir]
  snd_esti_ctype_labels = ctypes[esti_ctype_sec]

  # ARI(fst_esti_ctype_labels, true_ctype_labels)
  # table(fst_esti_ctype_labels, true_ctype_labels)

  ari = ARI(fst_esti_ctype_labels[true_ctype_labels_h != 'DBL'], true_ctype_labels_h[true_ctype_labels_h != 'DBL'])
  print(ari)
  # table(fst_esti_ctype_labels[true_ctype_labels_h != 'DBL'], true_ctype_labels_h[true_ctype_labels_h != 'DBL'])
  print('round')
  print(i)
}


table(fst_esti_ctype_labels[true_ctype_labels_h != 'DBL'], true_ctype_labels_h[true_ctype_labels_h != 'DBL'])
tb = table(fst_esti_ctype_labels[true_ctype_labels_h != 'DBL'], true_ctype_labels_h[true_ctype_labels_h != 'DBL'])

write.csv(tb, 'confusion_cell_line_high_PACS.csv')


##################################
## stratified sampling end
##################################




##################################
## stratified sampling -- low to high
##################################


  r_by_ct_low = get_r_by_ct_mat_pq(cell_type_set = cell_type_set,
                                   r_by_c = p_by_c,
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
                                       in_r_by_c = p_by_c_h[pk_sel,] )
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
p_by_c = p_by_c[diff_prop >= 0.05 & rowsd >= 0.015 ,] ## 188055 features left
p_by_c_h = p_by_c_h[diff_prop >= 0.05 & rowsd >= 0.015 ,] ## 188055 features left


# p_ct <- get_p_by_ct_mat(p_by_c = p_by_c,group_set = cell_type_set,group_lables = cell_type_labels,
#                         capturing_rate =  q_vec,max_p = 0.999, min_p = 0.001,adding_doublet = T)

# saveRDS(p_ct, 'p_by_ct_with_doublet.rds')
#
# p_ct <- readRDS('p_by_ct_with_doublet.rds')


#
# est_m <- estimate_type(p_by_t = p_ct, in_p_by_c = p_by_c, capturing_rate = q_vec, alpha = 1 )
#
# saveRDS(est_m, 'est_m_ArchR_R4.1.rds')




## test alternative doublet models
# p_ct_ave <- get_p_by_ct_mat_ave(p_by_c = p_by_c,group_set = cell_type_set,group_lables = cell_type_labels,
#                                 capturing_rate =  q_vec,max_p = 0.999, min_p = 0.001,adding_doublet = T)


p_ct_ave <- get_p_by_ct_mat_ave(p_by_c = p_by_c,group_set = cell_type_set,group_lables = cell_type_labels,
                                capturing_rate =  q_vec,max_p = 0.999, min_p = 0.001,adding_doublet = T)

save.image(file = "ArchR_whole.RData")

# saveRDS(p_ct_ave, 'p_ct_ave_new_filt.rds')

# p_by_c_dou <- p_by_c[,cell_type_labels == 'DBL']
# q_vec_dou <- q_vec[cell_type_labels == 'DBL']

# est_m_doub <- estimate_type(p_by_t = p_ct_ave, in_p_by_c = p_by_c_dou,
#                             capturing_rate = q_vec_dou, alpha = 1 )
#
# saveRDS(est_m_doub, 'est_m_ArchR_doub.rds')

# est_m_doub <- estimate_type(p_by_t = p_ct_ave, in_p_by_c = p_by_c,
#                             capturing_rate = q_vec, alpha = 1 )
est_m_doub <- estimate_type(p_by_t = p_ct_half, in_p_by_c = p_by_c,
                            capturing_rate = q_vec, alpha = 1 )


saveRDS(est_m_doub, 'est_m_ArchR_doub_half.rds')






