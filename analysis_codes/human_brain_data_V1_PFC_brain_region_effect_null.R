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
library('PICsnATAC')
library('PACS')



brain = readRDS('brain_mat_PFC_V1.rds')  # 434353  21781
meta_use = readRDS('meta_use_PFC_V1.rds')

## get dummy variable -- location effect, but keep all

xdummy <- meta_use %>% mutate(dummy=1) %>%
  spread(key=CellType,value=dummy, fill=0) %>% select(-RG)
xdummy <- xdummy %>% mutate(dummy=1) %>%
  spread(key=specimen,value=dummy, fill=0) %>% select(-GW20)
xdummy <- xdummy %>% mutate(dummy=1) %>%
  spread(key=area,value=dummy, fill=0) %>% select(-PFC)
xdummy <- xdummy %>% select(-uniqueID)
xdummy <- as.matrix(xdummy)
xdummy <- cbind(1,xdummy)
colnames(xdummy)[1] <- 'intercept'

## notice, the location effect is the last column


## get cell type labels
cell_type_labels <- meta_use$CellType
cell_type_set <- unique(cell_type_labels)
rm(meta_use)
gc()

r_by_ct = readRDS('brain_r_by_ct_PFC_V1.rds')

par_initial_one = rep(0.02, length = dim(xdummy)[2])
par_initial_one[1] = -0.05
par_initial_one[8] = 0

brain = brain[(1+(iii-1)*8000):(iii*8000),]
brain = as.matrix(brain)
gc()


no_area_para = estimate_parameters_null(r_by_c = brain,
               design_mat = xdummy,
               par_initial = par_initial_one,
               hold_zero = 8,
               cap_rate_vec = r_by_ct$q_vec_new,
               mc.cores = 1)

saveRDS(no_area_para, paste0('no_area_para_', iii,'.rds'))


