#BiocManager::install("TCGAbiolinks")
#BiocManager::install("EDASeq")

#Uncomment the following section for data preprocessing
# library(TCGAbiolinks)
# library(dplyr)
# library(DT)
# library(EDASeq)
# library(mice)

library(MLmetrics)
library(randomForest)
library(glmnet)
library(rpart)
library(nnet)
library(nnls)
library(foreach)
library(parallel)
library(doParallel)
library(genefilter)
library(abind)
library(locfit)
library(cluster)
library(mvtnorm)
library(mclust)
library(vscc)
library(openxlsx)


######### Molecular subtypes
#### Uncomment the following section for data preprocessing

# query <- GDCquery(project = "TCGA-LGG",
#                   data.category = "Gene expression",
#                   data.type = "Gene expression quantification",
#                   platform = "Illumina HiSeq",
#                   file.type  = "results",
#                   sample.type = c("Primary Tumor"),
#                   legacy = TRUE)
# GDCdownload(query)
# query$results[[1]]
#
#
#
# lgg.exp <- GDCprepare(query, save = TRUE, summarizedExperiment = TRUE, save.filename = "LGGIllumina_HiSeq.rda")
# genedata <- t(scale(assays(lgg.exp)[[1]]))
#
# gd <- order(apply(genedata, 1, var), decreasing=TRUE)[1:100]
# genedata <- genedata[, gd]
# write.table(genedata, paste0("~/Desktop/lgg_genedata.csv"), sep = ",")
#
# lgg.gbm.subtype <- TCGAquery_subtype(tumor = "lgg")


#For the actual lgg subtype
# lgg <- as.data.frame(lgg.gbm.subtype)
# lgg_num <- as.data.frame(sapply(1:ncol(lgg),  function(i) as.numeric(lgg[,i])))
# colnames(lgg_num) = colnames(lgg)

#outcome: number of months for survival

# outcome <- lgg_num$Mutation.Count
# test_indices <- which(is.na(outcome))
# y <- na.omit(outcome)
#
# full_data <- lgg_num[- test_indices, -20]

######

#miceMod <- mice(full_data, method="rf")
#imputed_data <- complete(miceMod)

#for multifaceted:
#data_full <- imputed_data[, -1]

#for (j in 1:ncol(data_full)){
#  data_full[,j] <- as.numeric(data_full[,j])
#}


#for lgg
# data_full <- imputed_data[, names(full_data) %in% c("Age..years.at.diagnosis.", "Karnofsky.Performance.Score",
#                                                     "Survival..months.","Percent.aneuploidy", "TERT.expression..log2.",
#                                                     "ABSOLUTE.purity", "ABSOLUTE.ploidy", "ESTIMATE.stromal.score",
#                                                     "ESTIMATE.immune.score", "ESTIMATE.combined.score", "Telomere.length.estimate.in.blood.normal..Kb.",
#                                                     "Telomere.length.estimate.in.tumor..Kb.")]
#
#
# data_full <- cbind(y, data_full)
#
# rest_full <- imputed_data[, !names(full_data) %in% c("Age..years.at.diagnosis.", "Karnofsky.Performance.Score",
#                                                     "Survival..months.","Percent.aneuploidy", "TERT.expression..log2.",
#                                                     "ABSOLUTE.purity", "ABSOLUTE.ploidy", "ESTIMATE.stromal.score",
#                                                     "ESTIMATE.immune.score", "ESTIMATE.combined.score", "Telomere.length.estimate.in.blood.normal..Kb.",
#                                                     "Telomere.length.estimate.in.tumor..Kb.")][,-1]
#


#factor_ind <- which(sapply(1:ncol(data_full), function(i) length(unique(data_full[,i]))) < 50)

#for (i in factor_ind){
#  data_full[, i]  <- as.factor(data_full[,i])
#}

#takeout <- which(sapply(1:ncol(data_full), function(i) sum(is.na(data_full[,i]))) > 0 )

#only if takeout > 0
#data_full <- data_full[, -takeout]

#write.table(data_full, paste0("~/Desktop/lgg_data_mutation.csv"), sep = ",")
#write.table(data_full, paste0("~/Desktop/lgg_rest_mutation.csv"), sep = ",")
#write.table(cbind(y, imputed_data), paste0("~/Desktop/lgg_imputed_mutation.csv"), sep = ",")

test_indices <- c(107, 154, 339)

ind_merged <- read.csv("/home/mr356/cluster_project/data_examples/lgg/ind_merged.csv")
ind_test <- read.csv("/home/mr356/cluster_project/data_examples/lgg/ind_test.csv")

data_full <- as.data.frame(read.csv("/home/mr356/cluster_project/data_examples/lgg_mutation/lgg_data_mutation.csv"))
genedata <- as.data.frame(read.csv("/home/mr356/cluster_project/data_examples/lgg_gene/lgg_genedata.csv"))

imputed_data <- as.data.frame(read.csv("/home/mr356/cluster_project/data_examples/lgg_mutation_fullmerged/lgg_imputed_mutation.csv"))
imputed_data <- imputed_data[, !(names(imputed_data) %in% c( "Pan.Glioma.RNA.Expression.Cluster", "Random.Forest.Sturm.Cluster",          
                                                             "Pan.Glioma.DNA.Methylation.Cluster", "IDH.specific.DNA.Methylation.Cluster", "Supervised.DNA.Methylation.Cluster" ))]


dfull <- cbind(imputed_data, genedata[-test_indices, ])

dfull <- cbind(imputed_data[,names(imputed_data) %in% c("Grade", "IDH.specific.RNA.Expression.Cluster")], genedata[-test_indices, ])
dfull[,1] <- as.factor(dfull[,1])

ind_merged <- c(ind_merged$ind_merged)
ind_test <- c(ind_test$ind_test)


randomforestfit <- function(data, ...){
  rf <- randomForest::randomForest(y ~ ., data = as.data.frame(data), ntree = 100, importance = TRUE, ...)
  rf
}

randomforestpredict <- function(data, newdata){
  as.vector(predict(data, newdata=as.data.frame(newdata)))
}

rcomb <- function(...) {
  args <- list(...)
  lapply(seq_along(args[[1]]), function(i)
    do.call('rbind', lapply(args, function(a) a[[i]])))
}

LogLoss <- function (y_pred, y_true){
  eps <- 1e-15
  y_pred <- pmax(pmin(y_pred, 1 - eps), eps)
  LogLoss <- -mean(y_true * log(y_pred) + (1 - y_true) * log(1 - y_pred))
  return(LogLoss)
}

absnorm <- function(vec, max.norm = FALSE){
  sgn <- sign(vec)
  vec <- abs(vec)
  if(max.norm){
    mvec <- max(vec)	
  } else {
    mvec <- min(vec)
  }
  
  vec <- abs(vec - mvec)
  sgn*(vec/sum(vec))
}

#clustering using mclust
create_clusters <- function(orig_data, X, ntest, k){
  change_order <- sample(nrow(X))
  
  X <- X[change_order, ]
  merged <- X[-c(1:ntest), ]
  #merged <- merged[sample(nrow(merged)), ]
  #cluster without using y
  
  #mclust
  vm <- vscc(x = merged[,-1], G = c(5, 9, 13), forcereduction = TRUE)
  print(vm$bestmodel$G)
  k2 <- vm$bestmodel$classification
  
  #k_means
  # k2 <- kmeans(merged[,-1], centers = k, nstart = 25, iter.max = 25)
  # if (k2$ifault==4) { 
  #   k2 = kmeans(merged[,-1], k2$centers, nstart = 25, iter.max = 25, algorithm="MacQueen")$cluster 
  # }else{
  #   k2 <- k2$cluster
  # }
  # 
  orig_data <- orig_data[change_order, ]
  training_data <- orig_data[-c(1:ntest), ]
  
  clusters_list <- lapply(split(seq_along(k2), k2), 
                          function(m, ind) m[ind,], m = training_data)[order(unique(k2))]
  len_clust <- sapply(clusters_list, function(i) nrow(i))
  wlc <- which(len_clust <= 2)
  if (any(wlc)){
    clusters_list <- clusters_list[-wlc]
  }
  nclusters <- length(clusters_list)
  #add test set
  clusters_list[[nclusters + 1]] <- orig_data[c(1:ntest), ]
  for (j in 1:length(clusters_list)){
    clusters_list[[j]][,1] <- as.factor(clusters_list[[j]][,1])
  }
  
  return(list(clusters_list = clusters_list))
}



#Using the true clusters from the paper
create_clusters_comb <- function(data_full, nsep, ind_merged, ind_test){
  ntest = 1
  ndat = nsep + ntest
  edat <-  vector("list", ndat)
  d <- (1:nrow(data_full))[sample(1:nrow(data_full))]
  x <- seq_along(d)
  d1 <- split(d, ceiling(x/(nrow(data_full)/ndat)))
  
  
  ind_merged <- do.call(c, d1[1:(length(d1) - 1)])
  ind_test = d1[[length(d1)]]
  
  merged <- data_full[ind_merged,]
  test <- data_full[ind_test, ]
  
  merged <- dfull[ind_merged,]
  test <- dfull[ind_test,]
  k2 <- merged$IDH.specific.RNA.Expression.Cluster
  
  clusters_list <- lapply(split(seq_along(k2), k2), #split indices by a
                          function(m, ind) m[ind, ], m = merged)[order(unique(k2))] #merged_choose
  len_clust <- sapply(clusters_list, function(i) nrow(i))
  wlc <- which(len_clust <= 5)
  if (any(wlc)){
    clusters_list <- clusters_list[-wlc]
  }
  nclusters <- length(clusters_list)
  
  for (t in 1:ntest){
    clusters_list[[(nclusters + t)]] <- test #test_choose
  }
  for (j in 1:length(clusters_list)){
    sc <- clusters_list[[j]]
    sc[is.na(sc)] <- 0
    clusters_list[[j]] <- as.data.frame(sc)
    colnames(clusters_list[[j]]) <- c("y", paste0("V", 1:(ncol(sc)-1)))
    
  }
  unique_ls <- lapply(1:length(clusters_list), function(i) unique(clusters_list[[i]]$y))
  for (t in  1:length(unique_ls)){
    if (length(unique_ls[[t]]) == 1){
      clusters_list  <- clusters_list[-t]
    }
  }
  return(list(clusters_list = clusters_list, ind_merged = ind_merged, ind_test = ind_test)) 
}

full_merged <- function(data_full, ind_merged, ind_test, modpred){
  ntest = 1 
  data <- data_full[ind_merged, ]
  data <- sapply(1:ncol(data), function(i) as.numeric(data[,i]))
  nvar <- ncol(data) - 1
  colnames(data) <- c("y", paste0("V", 1:nvar))
  
  test_set <- data_full[ind_test,]
  test_set <- sapply(1:ncol(test_set), function(i) as.numeric(test_set[,i]))
  colnames(test_set) <- c("y", paste0("V", 1:nvar))
  
  mod0 <- randomForest::randomForest(y ~ ., data , ntree = 500, importance = TRUE)
  mod0_importance <- order(importance(mod0)[,1], decreasing =  TRUE)
  full_preds <- modpred(mod0, newdata = test_set)
  cury <- as.numeric(as.character(test_set[,1]))
  rmse <- sqrt(mean((cury - full_preds)^2))
  return(list(rmse  = rmse, mod0_importance  = mod0_importance))
}

clusters_fit <- function(modfit, modpred, sd){
  edat = sd$clusters_list
  ind_merged = sd$ind_merged
  ind_test = sd$ind_test
  cols_4 = sd$cols_4
  
  ntest = 1
  
  ntrain = length(edat) - ntest
  
  mods <- vector("list", ntrain)
  mses <- matrix(NA, ntrain, ntrain)
  
  allpreds <- rf_importances <- vector("list", ntrain)
  
  #learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain])
  
  mod0 <- randomForest::randomForest(y ~ ., data = as.data.frame(matstack[sample(nrow(matstack)), ]), ntree = 500, importance = TRUE)
  
  mod0_importance <- order(importance(mod0)[,1], decreasing =  TRUE)
  
  for (j in 1:ntrain){
    mods[[j]] <- modfit(edat[[j]])
    rf_importances[[j]] <- order(importance(mods[[j]])[,1], decreasing =  TRUE)
    
    preds <- lapply(edat[1:ntrain], function(x){
      modpred(mods[[j]], newdata = x[, -1]) 
    })
    mses[j,] <- unlist(lapply(edat[1:ntrain], function(x){#cross validation within the training set
      newdata = x[, -1]
      preds <- modpred(mods[[j]], newdata = newdata) 
      mean((preds - x[,"y"])^2)}
    ))
    
    curpreds <- lapply(preds, as.numeric)
    allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
  }
  diag(mses) <- NA
  rfi <-  do.call(rbind, rf_importances[1:ntrain])
  
  matstack <- as.data.frame(matstack)
  # CS Weights
  tt <- apply(mses, 1, mean, na.rm = T) #removing the diagonal elements, takes the mean of each row 
  weights <- absnorm(sqrt(tt), max.norm = TRUE)
  nk <- unlist(lapply(edat, nrow)) #vector of number of rows in each dataset
  nwts <- absnorm(nk[1:ntrain])
  
  # Regression: stacked (intercept and no intercept)
  predstack <- do.call(rbind, allpreds)
  coefs_stack_noint <- nnls::nnls(predstack, as.numeric(as.character(matstack$y)))$x
  coefs_stack_int <- nnls::nnls(cbind(rep(1,nrow(predstack)),predstack), as.numeric(as.character(matstack$y)))$x
  coefs_stack_lasso <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 1, lower.limits = 0, intercept = T)))
  coefs_stack_ridge <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 0, lower.limits = 0, intercept = T)))
  #Just a safeguard against full collinearity, although I think we are OK with nnls now
  coefs_stack_noint[which(is.na(coefs_stack_noint))] <- 0
  coefs_stack_int[which(is.na(coefs_stack_int))] <- 0
  coefs_stack_lasso[which(is.na(coefs_stack_lasso))] <- 0
  coefs_stack_ridge[which(is.na(coefs_stack_ridge))] <- 0
  
  coefs_stack_noint_norm <- absnorm(coefs_stack_noint)
  
  # Regression: study-specific (intercept and no intercept)
  coefs_ss_noint <- mapply(function(x,y){nnls::nnls(y,as.numeric(as.character(x[,"y"])))$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_noint <- colMeans(do.call(rbind, coefs_ss_noint), na.rm = T)
  coefs_ss_int <- mapply(function(x,y){nnls::nnls(cbind(rep(1,nrow(y)),y),as.numeric(as.character(x[,"y"])))$x}, edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_int <- colMeans(do.call(rbind, coefs_ss_int), na.rm = T)
  coefs_ss_noint[which(is.na(coefs_ss_noint))] <- 0
  coefs_ss_int[which(is.na(coefs_ss_int))] <- 0
  coefs_ss_noint_norm <- absnorm(coefs_ss_noint)
  
  outmat <- matrix(NA, ntest, 12)
  colnames(outmat) <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                        "Stack_noint", "Stack_noint_norm", "Stack_int",
                        "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "Stack_ridge")
  
  for(i in (ntrain + 1):(length(edat))){
    merged <- modpred(mod0, newdata = edat[[i]][,-1])
    merged <- as.vector(sapply(merged, as.numeric))
    allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1]))
    allmod <- apply(allmod, 2, as.numeric)
    # Unweighted average
    unweighted <- colMeans(allmod)
    
    # sample size weighted
    sample_wtd <- apply(allmod, 2, function(x){sum(nwts*x)})
    
    # cross-study weighted 
    cs_wtd <- apply(allmod, 2, function(x){sum(weights*x)})
    
    # regression: stacked (noint, int, each normed) + lasso
    stack_noint <- apply(allmod, 2, function(x){sum(coefs_stack_noint*x)})
    stack_noint_norm <- apply(allmod, 2, function(x){sum(coefs_stack_noint_norm*x)})
    stack_int <- apply(allmod, 2, function(x){coefs_stack_int[1] + sum(coefs_stack_int[-1]*x)})
    stack_lasso <- apply(allmod, 2, function(x){coefs_stack_lasso[1] + sum(coefs_stack_lasso[-1]*x)})
    stack_ridge <- apply(allmod, 2, function(x){coefs_stack_ridge[1] + sum(coefs_stack_ridge[-1]*x)})
    
    # regression: study_specific (noint, int, noint normed) + lasso
    ss_noint <- apply(allmod, 2, function(x){sum(coefs_ss_noint*x)})
    ss_noint_norm <- apply(allmod, 2, function(x){sum(coefs_ss_noint_norm*x)})
    ss_int <- apply(allmod, 2, function(x){coefs_ss_int[1] + sum(coefs_ss_int[-1]*x)})
    
    cury <- as.numeric(as.character(as.data.frame(edat[[i]])[,"y"]))
    
    outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2), mean((cury - unweighted)^2), mean((cury - sample_wtd)^2), mean((cury - cs_wtd)^2),
                                  mean((cury - stack_noint)^2), mean((cury - stack_noint_norm)^2), mean((cury - stack_int)^2),
                                  mean((cury - ss_noint)^2), mean((cury - ss_noint_norm)^2), 
                                  mean((cury - ss_int)^2), mean((cury - stack_lasso)^2), 
                                  mean((cury - stack_ridge)^2)))
  }
  outmat <- (outmat - outmat[,1])/outmat[,1]*100
  return(list(outmat = colMeans(outmat), mod0_importance = mod0_importance))
}


rep.clusters_fit <- function(reps, modfit, modpred, ind_merged, ind_test){
  #####
  logfile <- paste0("outputFile","nonsparsity",".txt")
  writeLines(c(""), file(logfile,'w'))
  
  num.threads <- round(reps/5)# round up
  threads <- makeCluster(num.threads, outfile=logfile, setup_timeout = 0.5)
  registerDoParallel(threads)
  
  getDoParWorkers()
  #######
  
  results <- foreach(i = 1:reps, .combine = 'rcomb', .multicombine = TRUE, .export = ls(globalenv())) %dopar% {
    library(locfit)
    library(ClusterR)
    library(cluster)
    library(mvtnorm)
    library(MLmetrics)
    library(randomForest)
    library(glmnet)
    library(nnet)
    library(nnls)
    library(genefilter)
    library(mclust)
    library(vscc)
    
    #for cluster: 
    sd <- create_clusters_comb(data_full, nsep = 9, ind_merged, ind_test)
    cf <- clusters_fit(modfit, modpred, sd)
    errors_cluster <- cf$outmat
    #print(errors_cluster)
    chosen_cols <- c(sd$chosen_cols, rep(NA, 12 - length(sd$chosen_cols)))
    print(chosen_cols)
    ind_merged <- sd$ind_merged
    ind_test <- sd$ind_test
    errors_full <- full_merged(imputed_data, ind_merged, ind_test, modpred)
    errors_total <- c(errors_full$rmse, errors_cluster)
    importances <- cf$mod0_importance
    print(importances)
    errors_total <- (errors_total - errors_total[2])/errors_total[2]*100
    print(errors_total)
    return(list(errors_total = errors_total, chosen_cols  = chosen_cols, importances = importances))
  }
  closeAllConnections()
  errors_total = as.data.frame(results[[1]])
  chosen_cols = results[[2]]
  importances  =  results[[3]]
  
  colnames_total <- c("Full","Merged", "Unweighted", "Sample_Weighted","CS_Weighted", "Stack_noint", "Stack_noint_norm", "Stack_int",
                      "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "Stack_ridge")
  
  colnames(errors_total) <- colnames_total
  
  write.table(as.data.frame(errors_total), paste0("errors_total_50.csv"), sep = ",", col.names = colnames_total)
  write.table(chosen_cols, paste0("chosen_cols.csv"), sep = ",")
  write.table(importances, paste0("importances.csv"), sep = ",")
  
  
  return(list(errors_total = errors_total, chosen_cols = chosen_cols,  importances =  importances))   
}

rcf <- rep.clusters_fit(reps = 250, modfit = randomforestfit, modpred = randomforestpredict, ind_merged, ind_test)

