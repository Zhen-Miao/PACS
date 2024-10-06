## jobindex
jobindex <- as.numeric(Sys.getenv("LSB_JOBINDEX"))
if(is.na(jobindex)| jobindex ==0){jobindex = 1}

iii = jobindex

## load libraries
library('SummarizedExperiment')
library('dplyr')
library('Matrix')
library('tidyr')
library('parallel')
library('tictoc')
library('PACS')

figr = readRDS('pma_only_data_mat.rds')
meta_use = readRDS('pma_only_cell_meta.rds')

meta_use = meta_use[colnames(figr),]

meta_use = meta_use[,c('Donor','Condition','pairedLabel2')]
meta_use$Condition[meta_use$Condition == 'Control_1h'] = 0
meta_use$Condition[meta_use$Condition == 'PMA_1h'] = 1
meta_use$Condition[meta_use$Condition == 'PMA_6h'] = 2
meta_use$id = c(1:nrow(meta_use))

## filter cell types with too few cells
ct_sel = names(table(meta_use$pairedLabel2))[table(meta_use$pairedLabel2) >300]
cell_sel = rownames(meta_use)[meta_use$pairedLabel2 %in% ct_sel]

meta_use = meta_use[cell_sel,]
figr = figr[,cell_sel]


ctype_sel = names(table(meta_use$pairedLabel2))[iii]
cell_sel2 = rownames(meta_use)[meta_use$pairedLabel2 == ctype_sel]
meta_use = meta_use[cell_sel2,c(1,2,4)]
figr = figr[,cell_sel2]


## get dummy variable
xdummy <- meta_use %>% mutate(dummy=1) %>%
  spread(key=Donor,value=dummy, fill=0) %>% select(-Donor1)
xdummy <- xdummy %>% select(-id)

xdummy <- xdummy %>% mutate(time = Condition)
xdummy <- xdummy %>% select(-Condition)
xdummy$time = as.integer(xdummy$time)

xdummy <- as.matrix(xdummy)
xdummy <- cbind(1,xdummy)
colnames(xdummy)[1] <- 'intercept'



r_by_ct = readRDS('figr_r_by_ct_PMA.rds')
q_vec_new = r_by_ct$q_vec_new[colnames(figr)]



para_full = readRDS(file = paste('full_para_figr_PMA_',ctype_sel,'_cell.rds'))
para_no_area = readRDS(file = paste('partial_para_figr_PMA_',ctype_sel,'_cell.rds'))


conver_indi <- para_full[nrow(para_full),] == 1 &
  para_no_area[nrow(para_no_area),] == 1

figr = figr[conver_indi,]

figr_t = t(figr)
rm(figr)
gc()

para_filt = para_full[1:(nrow(para_full)-1),conver_indi]
para_no_area_filt = para_no_area[1:(nrow(para_no_area)-1),conver_indi]



q_vec = q_vec_new

compute_p_value = function(para_filt, para_no_area_filt,
                           figr_t, x_full,x_null, q_vec, df_test){
  lrt_area = vector(length = dim(figr_t)[2])
  names(lrt_area) = colnames(figr_t)
  n_features_iter = ceiling(ncol(figr_t) / 5000)

  for(i in 1:n_features_iter){
    from_i = 1 + (i-1) * 5000
    to_i = min(i * 5000, dim(figr_t)[2])

    c_by_r = as.matrix(figr_t[,from_i:to_i])
    features = colnames(c_by_r)

    lrt_area[features]  = compare_models(
      x_full = x_full,
      theta_estimated_full = para_filt[,from_i:to_i],
      x_null = x_null,
      theta_estimated_null = para_no_area_filt[,from_i:to_i],
      q_vec = q_vec,
      c_by_r = c_by_r,
      df_test = df_test,
      mc.cores = 1)
    print(i)
    rm(c_by_r)
    gc()
  }
  return(lrt_area)
}


compute_acc_change = function(kidney,  q_vec,cell_types,cell_type_set){

  n_features_iter = ceiling(nrow(kidney) / 5000)
  if(length(cell_type_set) == 3){
    fc_area = matrix(nrow = dim(kidney)[1], ncol = length(cell_type_set) - 1)
    rownames(fc_area) = rownames(kidney)
    colnames(fc_area) = cell_type_set[2:length(cell_type_set)]

    fc_area_no_missing = fc_area

    for(i in 1:n_features_iter){
      from_i = 1 + (i-1) * 5000
      to_i = min(i * 5000, dim(kidney)[1])
      n_cells = to_i - from_i + 1

      c_by_r = as.matrix(kidney[from_i:to_i,])
      features = rownames(c_by_r)
      cap_mat = 1 / q_vec[cell_types == cell_type_set[1]]
      cap_mat = matrix(rep(cap_mat, n_cells), nrow = n_cells, byrow = T)

      cap_mat_2 = 1 / q_vec[cell_types == cell_type_set[2]]
      cap_mat_2 = matrix(rep(cap_mat_2, n_cells), nrow = n_cells, byrow = T)

      cap_mat_3 = 1 / q_vec[cell_types == cell_type_set[3]]
      cap_mat_3 = matrix(rep(cap_mat_3, n_cells), nrow = n_cells, byrow = T)

      rm1 <- rowSums(c_by_r[,cell_types == cell_type_set[1] ] * cap_mat) /
        sum(cell_types == cell_type_set[1])
      rm2 <- rowSums(c_by_r[,cell_types == cell_type_set[2] ] * cap_mat_2) /
        sum(cell_types == cell_type_set[2])
      rm3 <- rowSums(c_by_r[,cell_types == cell_type_set[3] ] * cap_mat_3) /
        sum(cell_types == cell_type_set[3])

      ## sequencing depth corrected log fold change
      fc_area[features,1] = rm2 - rm1
      fc_area[features,2] = rm3 - rm1

      rm1_unc = rowMeans(c_by_r[,cell_types == cell_type_set[1]])
      rm2_unc = rowMeans(c_by_r[,cell_types == cell_type_set[2]])
      rm3_unc = rowMeans(c_by_r[,cell_types == cell_type_set[3]])

      fc_area_no_missing[features,1] = rm2_unc -  rm1_unc
      fc_area_no_missing[features,2] = rm3_unc -  rm1_unc

      print(i)
      rm(c_by_r)
      gc()
    }
    return(list(acc_change = fc_area, acc_change_no_missing = fc_area_no_missing ))
  }else{
    fc_area = vector(length = dim(kidney)[1])
    names(fc_area) = rownames(kidney)

    fc_area_no_missing = fc_area

    for(i in 1:n_features_iter){
      from_i = 1 + (i-1) * 5000
      to_i = min(i * 5000, dim(kidney)[1])
      n_cells = to_i - from_i + 1

      c_by_r = as.matrix(kidney[from_i:to_i,])
      features = rownames(c_by_r)

      cap_mat = 1 / q_vec[cell_types == cell_type_set[1]]
      cap_mat = matrix(rep(cap_mat, n_cells), nrow = n_cells, byrow = T)

      cap_mat_2 = 1 / q_vec[cell_types == cell_type_set[2]]
      cap_mat_2 = matrix(rep(cap_mat_2, n_cells), nrow = n_cells, byrow = T)

      rm1 <- rowSums(c_by_r[,cell_types == cell_type_set[1] ] * cap_mat) /
        sum(cell_types == cell_type_set[1])
      rm2 <- rowSums(c_by_r[,cell_types == cell_type_set[2] ] * cap_mat_2) /
        sum(cell_types == cell_type_set[2])

      ## sequencing depth corrected log fold change
      fc_area[features] = rm2 - rm1

      rm1_unc = rowMeans(c_by_r[,cell_types == cell_type_set[1]])
      rm2_unc = rowMeans(c_by_r[,cell_types == cell_type_set[2]])

      fc_area_no_missing[features] = rm2_unc -  rm1_unc

      print(i)
      rm(c_by_r)
      gc()
    }
    return(list(acc_change = fc_area, acc_change_no_missing = fc_area_no_missing ))
  }

}

fc = compute_acc_change(t(figr_t), q_vec, cell_types = meta_use$time,
                         cell_type_set = unique(meta_use$time))

saveRDS(fc, paste('acc_change_figr_PMA_',ctype_sel,'_cell.rds'))


lrt_area = compute_p_value(para_filt = para_filt,
                           para_no_area_filt = para_no_area_filt,
                           figr_t = figr_t, x_full = xdummy,
                           x_null = xdummy, q_vec, df_test = 1)


saveRDS(lrt_area, paste('lrt_time_p_valuesfigr_PMA_',ctype_sel,'_cell.rds'))

quantile(lrt_area)


