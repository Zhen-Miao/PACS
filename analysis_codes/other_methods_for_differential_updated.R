## all other methods for differential test
library(Rfast)
library('presto')
library('edgeR')

##################################################################
## The functions in this file are only for comparing with other methods
## we did not use functions here in the PACS workflow
##################################################################


###################################
## -- edgeR method
## default recommended by snapATAC
# 0.4 for human, 0.1 for mouse
#######################################
snapATAC_method <- function(data_matrix_pos, data_matrix_neg,bcv){
  if("dgCMatrix" %in% class(data_matrix_pos)){
    data_use = cbind(Matrix::rowSums(data_matrix_pos), Matrix::rowSums(data_matrix_neg))
  }else{
    data_use = cbind(Rfast::rowsums(data_matrix_pos), Rfast::rowsums(data_matrix_neg))
  }
  data_use = as.data.frame(data_use)
  group <- factor(c(1,2));
  design <- model.matrix(~group);
  y <- DGEList(counts=data_use, group=group);
  p_val_vec <- exactTest(y, dispersion=bcv^2)$table$PValue
  return(p_val_vec)
}

######################################################################
## the same thing, but accelerated
######################################################################
seurat_method2_subsample <- function(data_matrix_pos, data_matrix_neg, peak_region_fragments){
  data.use <- cbind(data_matrix_pos, data_matrix_neg)
  group.info <- c(rep(0,times = dim(data_matrix_pos)[2]),
                  rep(1,times = dim(data_matrix_neg)[2]) )

  # group.info <- as.factor(group.info)
  # peak_region_fragments <- Rfast::colsums(data.use)
  ## data.use should be feature by cell matrix
  n_features <- dim(data.use)[1]
  p_val <- vector(length = n_features)
  for(i in 1:n_features){
    X1 = cbind(data.use[i, ],peak_region_fragments)
    # model.data <- data.frame(GENE = data.use[i, ], group.info = group.info,
    #                          peak_region_fragments = peak_region_fragments)
    # model1 <- glm(formula = group.info ~ GENE + peak_region_fragments,
    #               data = model.data, family = binomial)
    # model2 <- glm(formula = group.info ~ peak_region_fragments,
    #               data = model.data, family = binomial)
    model1 <- glm_logistic(x = X1, y = group.info)
    model2 <- glm_logistic(x = peak_region_fragments,y = group.info)
    p_val[i] <- pchisq(model2$devi - model1$devi, df = 1,lower.tail = F)
  }
  return(p_val)
}


seurat_method3_subsample <- function(data_matrix, x_null, group.info){

  ## data_matrix should be feature by cell matrix
  n_features <- dim(data_matrix)[1]
  p_val <- vector(length = n_features)
  for(i in 1:n_features){
    X1 = cbind(x_null,data_matrix[i, ])

    model1 <- glm_logistic(x = X1, y = group.info)
    model2 <- glm_logistic(x = x_null,y = group.info)
    p_val[i] <- pchisq(model2$devi - model1$devi, df = 1,lower.tail = F)
  }
  return(p_val)
}

#############################################################
## this function is only for evaluation of the Seurat method
## for the test of its ability to detect significant effect
## of other covariates
#############################################################
seurat_method3_other_factors <- function(data_matrix, x_null, group_info, batch_info){

  ## data_matrix should be feature by cell matrix
  n_features <- dim(data_matrix)[1]
  p_val <- vector(length = n_features)
  for(i in 1:n_features){
    X_null_AC = cbind(x_null,data_matrix[i, ])
    X_1 = cbind(X_null_AC, batch_info)

    model1 <- glm_logistic(x = X_1, y = group_info)
    model2 <- glm_logistic(x = X_null_AC,y = group_info)
    p_val[i] <- pchisq(model2$devi - model1$devi, df = 1,lower.tail = F)
  }
  return(p_val)
}


# ## the same thing, but accelerated
# seurat_method2_subsample_change_response <- function(data_matrix_pos, data_matrix_neg, peak_region_fragments){
#   data.use <- cbind(data_matrix_pos, data_matrix_neg)
#   group.info <- c(rep('A',times = dim(data_matrix_pos)[2]),
#                   rep('B',times = dim(data_matrix_neg)[2]) )
#   group.info <- as.factor(group.info)
#   # peak_region_fragments <- Rfast::colsums(data.use)
#   ## data.use should be feature by cell matrix
#   n_features <- dim(data.use)[1]
#   p_val <- vector(length = n_features)
#   for(i in 1:n_features){
#     model.data <- data.frame(GENE = data.use[i, ], group.info = group.info,
#                              peak_region_fragments = peak_region_fragments)
#     model1 <- glm(formula = GENE ~ group.info + peak_region_fragments,
#                   data = model.data, family = binomial)
#     model2 <- glm(formula = GENE ~ peak_region_fragments,
#                   data = model.data, family = binomial)
#     p_val[i] <- pchisq(model2$deviance - model1$deviance, df = 1,lower.tail = F)
#   }
#   return(p_val)
# }

######################################################################
## fisher
######################################################################
fisher_method <- function(
  data_matrix_pos, data_matrix_neg
){
  n_cell_pos <- dim(data_matrix_pos)[2]
  n_cell_neg <- dim(data_matrix_neg)[2]
  n_features <- dim(data_matrix_pos)[1]
  p_val <- vector(length = n_features)


  if("dgCMatrix" %in% class(data_matrix_pos)){
    data_use1 = cbind(Matrix::rowSums(data_matrix_pos), Matrix::rowSums(data_matrix_neg))
  }else{
    data_use1 = cbind(Rfast::rowsums(data_matrix_pos), Rfast::rowsums(data_matrix_neg))
  }

  data_use2 = cbind(n_cell_pos-data_use1[,1], n_cell_neg - data_use1[,2])
  data_use <- cbind(data_use1, data_use2)

  for(i in 1:n_features){
    x_mat = matrix(data_use[i,], nrow = 2)
    p_val[i] <- fisher.test(x_mat)$p.value
  }

  return(p_val)

}

## archR functionality
"%ni%" <- function(x, table) !(match(x, table, nomatch = 0) > 0)


getQuantiles <- function(v = NULL, len = length(v)){
  if(length(v) < len){
    v2 <- rep(0, len)
    v2[seq_along(v)] <- v
  }else{
    v2 <- v
  }
  p <- trunc(rank(v2))/length(v2)
  if(length(v) < len){
    p <- p[seq_along(v)]
  }
  return(p)
}

computeKNN <- function(
  data = NULL,
  query = NULL,
  k = 50,
  includeSelf = FALSE,
  ...
){
  # .validInput(input = data, name = "data", valid = c("dataframe", "matrix"))
  # .validInput(input = query, name = "query", valid = c("dataframe", "matrix"))
  # .validInput(input = k, name = "k", valid = c("integer"))
  # .validInput(input = includeSelf, name = "includeSelf", valid = c("boolean"))
  if(is.null(query)){
    query <- data
    searchSelf <- TRUE
  }else{
    searchSelf <- FALSE
  }
  # .requirePackage("nabor", source = "cran")
  if(searchSelf & !includeSelf){
    knnIdx <- nabor::knn(data = data, query = query, k = k + 1, ...)$nn.idx
    knnIdx <- knnIdx[,-1,drop=FALSE]
  }else{
    knnIdx <- nabor::knn(data = data, query = query, k = k, ...)$nn.idx
  }
  knnIdx
}


matchBiasCellGroups <- function(
  input = NULL,
  groups = NULL,
  useGroups = NULL,
  bgdGroups = NULL,
  bias = NULL,
  k = 100,
  n = 500,
  bufferRatio = 0.8
){

  #Summary Function
  .summarizeColStats <- function(m = NULL, name = NULL){
    med <- apply(m, 2, median)
    mean <- colMeans(m)
    sd <- apply(m, 2, sd)
    loQ <- apply(m, 2, function(x) quantile(x, 0.25))
    hiQ <- apply(m, 2, function(x) quantile(x, 0.75))
    summaryDF <- t(data.frame(
      median = med,
      mean = mean,
      sd = sd,
      lowerQuartile = loQ,
      upperQuartile = hiQ
    )) %>% data.frame
    colnames(summaryDF) <- colnames(m)
    if(!is.null(name)){
      summaryDF$name <- name
    }
    summaryDF
  }



  #Make sure input is dataframe
  input <- data.frame(input)

  #Norm using input string ie log10(nfrags)
  inputNorm <- lapply(seq_along(bias), function(x){
    plyr::mutate(input, o=eval(parse(text=bias[x])))$o
  }) %>% Reduce("cbind", .) %>% data.frame
  rownames(inputNorm) <- rownames(input)

  #Quantile Normalization
  inputNormQ <- lapply(seq_len(ncol(inputNorm)), function(x){
    getQuantiles(inputNorm[,x])
  }) %>% Reduce("cbind", .) %>% data.frame
  rownames(inputNormQ) <- rownames(input)

  #Add Colnames
  colnames(inputNorm) <- bias
  colnames(inputNormQ) <- bias

  if(is.null(useGroups)){
    useGroups <- gtools::mixedsort(unique(paste0(groups)))
  }

  if(is.null(bgdGroups)){
    bgdGroups <- gtools::mixedsort(unique(paste0(groups)))
  }

  stopifnot(all(useGroups %in% unique(paste0(groups))))
  stopifnot(all(bgdGroups %in% unique(paste0(groups))))

  #Get proportion of each group
  prob <- table(groups) / length(groups)
  bgdProb <- prob[which(names(prob) %in% bgdGroups)] / sum(prob[which(names(prob) %in% bgdGroups)])

  #pb <- txtProgressBar(min=0,max=100,initial=0,style=3)
  matchList <- lapply(seq_along(useGroups), function(x){

    #setTxtProgressBar(pb,round(x*100/length(useGroups),0))

    #############
    # Organize
    #############
    groupx <- useGroups[x]
    idx <- which(names(bgdProb) == groupx)
    if(length(idx) > 0 & length(idx) != length(bgdProb)){
      bgdProbx <- bgdProb[-idx]/sum(bgdProb[-idx])
    }else{
      bgdProbx <- bgdProb
    }

    idF <- which(groups == groupx)

    if(all(length(idF) * bgdProbx < 1)){
      if(length(idF) < length(bgdProbx)){
        bgdProbx <- bgdProbx[sample(names(bgdProbx), floor(length(idF) * bufferRatio))]
        bgdProbx[1:length(bgdProbx)] <- rep(1/length(bgdProbx), length(bgdProbx))
      }
    }

    idB <- which(groups %in% names(bgdProbx))

    if(k > length(idB)){
      k2 <- length(idB)
    }else{
      k2 <- k
    }

    knnx <- computeKNN(inputNormQ[idB, ,drop=FALSE], inputNormQ[idF, ,drop=FALSE], k = k2)
    sx <- sample(seq_len(nrow(knnx)), nrow(knnx))

    minTotal <- min(n, length(sx) * bufferRatio)
    nx <- sort(floor(minTotal * bgdProbx))

    ###############
    # ID Matching
    ###############
    idX <- c()
    idY <- c()
    it <- 0

    if(any(nx <= 0)){
      nx[which(nx <= 0)] <- Inf
      nx <- sort(nx)
    }

    while(it < length(sx) & length(idX) < minTotal){

      it <- it + 1
      knnit <- knnx[sx[it],]
      groupit <- match(groups[idB][knnit],names(nx))
      selectUnique <- FALSE
      selectit <- 0
      oit <- order(groupit)

      while(!selectUnique){
        selectit <- selectit + 1
        itx <- which(oit==selectit)
        cellx <- knnit[itx]
        groupitx <- groupit[itx]
        if(is.infinite(nx[groupitx])){
          if(selectit == k2){
            itx <- NA
            cellx <- NA
            selectUnique <- TRUE
          }
        }else{
          if(cellx %ni% idY){
            selectUnique <- TRUE
          }
          if(selectit == k2){
            itx <- NA
            cellx <- NA
            selectUnique <- TRUE
          }
        }
      }

      if(!is.na(itx)){
        idX <- c(idX, sx[it])
        idY <- c(idY, cellx)
        nx[groupitx] <- nx[groupitx] - 1
        if(any(nx <= 0)){
          nx[which(nx <= 0)] <- Inf
          nx <- sort(nx)
        }
      }

      if(all(is.infinite(nx))){
        it <- length(sx)
      }

    }

    #####################
    # Convert Back to Normal Indexing
    #####################
    idX <- seq_len(nrow(inputNormQ))[idF][idX]
    idY <- seq_len(nrow(inputNormQ))[idB][idY]

    #####################
    # Matching Stats Groups
    #####################
    estbgd <- sort(floor(minTotal * bgdProbx))
    obsbgd <- rep(0, length(estbgd))
    names(obsbgd) <- names(estbgd)
    tabGroups <- table(groups[idY])
    obsbgd[names(tabGroups)] <- tabGroups
    estbgdP <- round(100 * estbgd / sum(estbgd),3)
    obsbgdP <- round(100 * obsbgd / sum(obsbgd),3)

    #####################
    # Matching Stats Bias Norm Values
    #####################
    forBias <- .summarizeColStats(inputNorm[idX,,drop=FALSE], name = "foreground")
    bgdBias <- .summarizeColStats(inputNorm[idY,,drop=FALSE], name = "background")

    out <- list(
      cells = idX,
      bgd = idY,
      summaryCells = forBias,
      summaryBgd = bgdBias,
      bgdGroups = rbind(estbgd, obsbgd),
      bgdGroupsProbs = rbind(estbgdP, obsbgdP),
      corbgdGroups = suppressWarnings(cor(estbgdP, obsbgdP)),
      n = length(sx),
      p = it / length(sx),
      group = groupx,
      k = k2
    )

    return(out)

  }) %>% SimpleList
  names(matchList) <- useGroups

  outList <- SimpleList(
    matchbgd = matchList,
    info = SimpleList(
      cells = rownames(input),
      groups = groups,
      biasNorm = inputNorm,
      biasNormQ = inputNormQ
    )
  )

  return(outList)

}


######################################################################
## archR
######################################################################
archR_method <- function(data_matrix_pos_count, data_matrix_neg_count){

  n_cells_pos <- dim(data_matrix_pos_count)[2]
  n_cells_neg <- dim(data_matrix_neg_count)[2]
  colnames(data_matrix_pos_count) <- as.character(c(1:n_cells_pos))
  colnames(data_matrix_neg_count) <- as.character(c(1:n_cells_neg))
  n_fragments <- c(Matrix::colSums(data_matrix_pos_count),
                   Matrix::colSums(data_matrix_neg_count))


  ## first, match the background cells with ArchR functions
  ## all input identical to ArchR default parameters
  colDat <- data.frame(logn_frag = log10(n_fragments))
  groups <- c( rep('A', times = n_cells_pos),
               rep('B', times = n_cells_neg))
  useGroups <- 'A'
  bgdGroups <- 'B'
  bias = "logn_frag"
  #####################
  ### check the default
  #####################
  k = 100
  bufferRatio = 0.8
  maxCells = 5000 ### this should be okay

  matchObj <- matchBiasCellGroups(
    input = colDat,
    groups = groups,
    useGroups = useGroups,
    bgdGroups = bgdGroups,
    bias = bias,
    k = k,
    n = maxCells,
    bufferRatio = bufferRatio
  )

  print('matchObj finished')

  matchx <- matchObj[[1]][[useGroups]]
  cellsx <- matchObj[[2]]$cells[matchx$cells]
  bgdx <- matchObj[[2]]$cells[matchx$bgd]


  ######## ------- show the length of matched cells ----
  print('length of matched cells')
  print(length(bgdx))
  print(head(bgdx))
  print(max(as.numeric(bgdx)))

  mat1 <- data_matrix_pos_count[,as.numeric(cellsx)]
  mat2 <- data_matrix_neg_count[,as.numeric(bgdx) - n_cells_pos]

  df <- wilcoxauc(cbind(mat1,mat2),
                  c(rep("Top", ncol(mat1)),
                    rep("Bot", ncol(mat2))))
  df <- df[which(df$group=="Top"),]
  return(df$pval)
}



archR_method_custom <- function(data_matrix_pos_count,
                                data_matrix_neg_count,
                                k = 20,
                                bufferRatio = 0.5,
                                maxCells = 5000){

  n_cells_pos <- dim(data_matrix_pos_count)[2]
  n_cells_neg <- dim(data_matrix_neg_count)[2]
  colnames(data_matrix_pos_count) <- as.character(c(1:n_cells_pos))
  colnames(data_matrix_neg_count) <- as.character(c(1:n_cells_neg))
  n_fragments <- c(Matrix::colSums(data_matrix_pos_count),
                   Matrix::colSums(data_matrix_neg_count))


  ## first, match the background cells with ArchR functions
  ## all input identical to ArchR default parameters
  colDat <- data.frame(logn_frag = log10(n_fragments))
  groups <- c( rep('A', times = n_cells_pos),
               rep('B', times = n_cells_neg))
  useGroups <- 'A'
  bgdGroups <- 'B'
  bias = "logn_frag"


  matchObj <- matchBiasCellGroups(
    input = colDat,
    groups = groups,
    useGroups = useGroups,
    bgdGroups = bgdGroups,
    bias = bias,
    k = k,
    n = maxCells,
    bufferRatio = bufferRatio
  )

  print('matchObj finished')

  matchx <- matchObj[[1]][[useGroups]]
  cellsx <- matchObj[[2]]$cells[matchx$cells]
  bgdx <- matchObj[[2]]$cells[matchx$bgd]


  ######## ------- show the length of matched cells ----
  print('length of matched cells')
  print(length(bgdx))
  print(head(bgdx))
  print(max(as.numeric(bgdx)))

  mat1 <- data_matrix_pos_count[,as.numeric(cellsx)]
  mat2 <- data_matrix_neg_count[,as.numeric(bgdx) - n_cells_pos]

  df <- wilcoxauc(cbind(mat1,mat2),
                  c(rep("Top", ncol(mat1)),
                    rep("Bot", ncol(mat2))))
  df <- df[which(df$group=="Top"),]
  return(df$pval)
}





