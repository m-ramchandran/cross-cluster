library(locfit)
library(fungible)
library(mvtnorm)
library(ClusterR)
library(cluster)
library(MLmetrics)
library(randomForest)
library(glmnet)
library(nnet)
library(nnls)
library(doParallel)
library(foreach)
library(plyr)
library(BiocGenerics)
library(clusterGeneration)


rcomb <- function(...) {
  args <- list(...)
  lapply(seq_along(args[[1]]), function(i)
    do.call('rbind', lapply(args, function(a) a[[i]])))
}

norm_vec <- function(x) sqrt(sum(x^2))

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

randomforestfit <- function(data, ...){
  rf <- randomForest::randomForest(y ~ ., data = as.data.frame(data), ntree = 100, importance = TRUE, ...)
  rf
}

randomforestpredict <- function(data, newdata){
  as.vector(predict(data, newdata=newdata))
}

#####
# Setup 1: Simulating data using the monte function - non-gaussian clusters

sim_data <- function(nstudies, ncoef, ntest){
  studies_list  <- vector("list", nstudies)
  nchoose <- 10
  #general predictor-outcome rule:
  coefs <- sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5)))
  vars <- sample(1:ncoef, nchoose)
  icoefs <- c(4, 1.8)
  
  #norm.betas <- norm_vec(c(coefs,icoefs))
  #coefs <- coefs/norm.betas
  #icoefs <- icoefs/norm.betas
  
  cormat <- matrix(.80, ncoef, ncoef)
  diag(cormat) <- rep(1, ncoef)
  in.cor.list <- replicate(nstudies, cormat, simplify=FALSE)
  
  
  m <- monte(seed = sample(1:1000, 1), nvar = ncoef, nclus = nstudies,
             clus.size = floor(runif(nstudies, (2000/(nstudies - ntest)), (2500/(nstudies - ntest)))), eta2 = runif(ncoef, .5, 1),
             cor.list = NULL, random.cor = FALSE, skew.list = NULL,
             kurt.list = NULL, secor = NULL, compactness = NULL,
             sortMeans = FALSE)
  
  for (i in 1:nstudies){
    curcoefs <- sapply(coefs, function(x){runif(1, x - .5, x + .5)})
    data_i <- scale(as.data.frame(m$data[m$data[,1] == i, ][,-1]))
    
    
    #generate outcome
    #scaled here, but the original variables are left unscaled for the clustering step
    y <- as.matrix((data_i[,vars]) %*% curcoefs) +
      icoefs[1]*(data_i[,vars[1]])^2 + icoefs[2]*(data_i[,vars[2]])^2 +
      cbind(rnorm(nrow(data_i))) # Added noise
    
    studies_list[[i]] <- as.data.frame(cbind(y, data_i))
    colnames(studies_list[[i]]) <- c("y", paste0("V", 1:ncoef))
  }
  return(list(studies_list = studies_list))
}

#####
# Setup 2: Simulating data using the gaussian mixture model paradigm

# sim_data <- function(nstudies, ncoef, ntest){
#   studies_list  <- vector("list", nstudies)
#   nchoose <- 10
#   #nchoose <- ncoef
#   #general predictor-outcome rule: 
#   coefs <- sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5))) 
#   vars <- sample(1:ncoef, nchoose)
#   icoefs <- c(4, 1.8)
#   n.noise <- 5
#   #If linear outcome model:
#   #icoefs <- c(0, 0)
#   
#   #norm.betas <- norm_vec(c(coefs,icoefs))
#   #coefs <- coefs/norm.betas
#   #cicoefs <- icoefs/norm.betas
#   
#   m <- genRandomClust(numClust = nstudies - ntest,
#                       sepVal=0.5,
#                       numNonNoisy=ncoef - n.noise,
#                       numNoisy=n.noise,
#                       numOutlier=100,
#                       numReplicate=1,
#                       fileName="test",
#                       clustszind=2,
#                       clustSizeEq=200,
#                       rangeN=c(300,350),
#                       clustSizes=NULL,
#                       covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),
#                       rangeVar=c(1, 10))
#   
#   m2 <- genRandomClust(numClust = 2,
#                        sepVal=0.8,
#                        numNonNoisy=ncoef - n.noise,
#                        numNoisy=n.noise,
#                        numOutlier=50,
#                        numReplicate=ntest,
#                        fileName="test",
#                        clustszind=2,
#                        clustSizeEq=50,
#                        rangeN=c(500,600),
#                        clustSizes=NULL,
#                        covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),
#                        rangeVar=c(1, 10))
#   
#   outliers <- scale(as.data.frame(m$datList)[which(m$memList[[1]] == 0), ])
#   outliers_list <- vector("list", nstudies)
#   
#   d <- 1:nrow(outliers)
#   x <- seq_along(d)
#   d1 <- split(d, ceiling(x/(nrow(outliers)/nstudies)))
#   for (i in 1:length(d1)){
#     outliers_list[[i]] <- as.data.frame(outliers[d1[[i]],])
#   }
#   
#   for (i in 1:nstudies){
#     curcoefs <- sapply(coefs, function(x){runif(1, x - .5, x + .5)})
#     
#     if (i <= (nstudies - ntest)){
#       data_i <- as.data.frame(m$datList)[which(m$memList[[1]] == i), ]
#       data_i  <- scale(rbind(data_i, outliers_list[[i]]))
#     }
#     else{
#       data_i <-scale(as.data.frame(m2$datList[[i - (nstudies - ntest)]]))
#     }
#     #data_i <- scale(as.data.frame(m$X[m$id == i, ]))
#     
#     
#     #generate outcome
#     #scaled here, but the original variables are left unscaled for the clustering step
#     y <- as.matrix((data_i[,vars]) %*% curcoefs) + 
#       icoefs[1]*(data_i[,vars[1]])^2 + icoefs[2]*(data_i[,vars[2]])^2 +
#       cbind(rnorm(nrow(data_i))) # Added noise
#     #new_y <- ifelse(y >= quantile(y, runif(1, .5, .9)), 1, 0)
#     
#     #studies_list[[i]] <- as.data.frame(cbind(as.factor(new_y), data_i))
#     #colnames(studies_list[[i]]) <- c("y", paste0("V", 1:ncoef))
#     
#     studies_list[[i]] <- as.data.frame(cbind(y, data_i))
#     colnames(studies_list[[i]]) <- c("y", paste0("V", 1:ncoef))
#   }
#   return(list(studies_list = studies_list))
# }

silhouette_score <- function(K_ls, merged){
  ss_ls <- c()
  km_ls <- vector("list", length(K_ls))
  for (i  in 1:length(K_ls)){
    k <- K_ls[i]
    km <- kmeans(merged, centers = k, nstart=25)
    ss <- mean(silhouette(km$cluster, dist(merged))[,3])
    ss_ls <- c(ss_ls, ss)
    km_ls[[i]] <- km$cluster
  }
  index_max <- which.max(ss_ls)
  #print(index_max)
  return(list(km_ind = km_ls[[index_max]], index_max = index_max))
}


create_clusters <- function(studies_list, ntest, ncoef){
  merged <- do.call(rbind, studies_list[1:(length(studies_list) - ntest)])
  merged <- merged[sample(nrow(merged)), ]
  #cluster without using y
  ss <- silhouette_score(2:(floor(nrow(merged)/30)), merged[,-1])
  #k2 <- kmeans(merged[,-1], centers = nclusters, nstart = 25)$cluster
  k2 <- ss$km_ind
  index_max <- ss$index_max
  clusters_list <- lapply(split(seq_along(k2), k2), #split indices by a
                          function(m, ind) m[ind,], m = merged)[order(unique(k2))]
  len_clust <- sapply(clusters_list, function(i) nrow(i))
  wlc <- which(len_clust <= 2)
  if (any(wlc)){
    clusters_list <- clusters_list[-wlc]
  }
  nclusters <- length(clusters_list)
  for (t in 1:ntest){
    clusters_list[[(nclusters + t)]] <- studies_list[[(length(studies_list) - ntest + t)]]
  }
  return(list(clusters_list = clusters_list, index_max = index_max))
}


create_random <- function(studies_list, ntest, ncoef, nsep){
  ntr <- length(studies_list) - ntest
  edat <-  vector("list", nsep + ntest)
  merged <- do.call(rbind, studies_list[1:ntr])
  merged <- merged[sample(nrow(merged)), ]
  
  
  d <- (1:nrow(merged))[sample(1:nrow(merged))]
  x <- seq_along(d)
  d1 <- split(d, ceiling(x/(nrow(merged)/nsep)))
  for (i in 1:length(d1)){
    edat[[i]] <- as.data.frame(merged[d1[[i]],])
  }
  for (t in 1:ntest){
    edat[[nsep + t]] <- as.data.frame(studies_list[[(length(studies_list) - ntest + t)]])
  }
  return(edat)
}

#cluster_ind = 1 if running the algorithm on the k-means clusters, 2 if randomly generating substudies, 3 if just running on the original set of simulated studies
clusters_fit <- function(modfit, modpred, ndat, ncoef, ntest, studies_list, cluster_ind, nsep){
  if (cluster_ind == 1){ #K-means clustering
    cc <- create_clusters(studies_list, ntest, ncoef)
    edat <- cc$clusters_list
    index_max <- cc$index_max
  }
  else if (cluster_ind == 2){ #random 'pseudo-studies'
    edat <- create_random(studies_list, ntest, ncoef, nsep)
  }
  else if (cluster_ind == 3){ #multistudy learning KNOWING the true clusters
    edat <- studies_list
  }
  
  ntrain = length(edat) - ntest
  
  mods <- vector("list", ntrain)
  mses <- matrix(NA, ntrain, ntrain)
  
  allpreds <- vector("list", ntrain)
  
  #learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain])
  
  #mod0 <- modfit(matstack[sample(nrow(matstack)), ])
  mod0 <- randomForest::randomForest(y ~ ., data = as.data.frame(matstack[sample(nrow(matstack)), ]), ntree = 500, importance = TRUE)
  
  
  
  for (j in 1:ntrain){
    mods[[j]] <- modfit(edat[[j]])
    
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
  coefs_ss_lasso <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,#cbind(rep(1,nrow(y)),y),
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 1, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_ridge <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,#cbind(rep(1,nrow(y)),y),
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 0, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_lasso <- colMeans(do.call(rbind, coefs_ss_lasso), na.rm = T)
  coefs_ss_ridge <- colMeans(do.call(rbind, coefs_ss_ridge), na.rm = T)
  coefs_ss_noint[which(is.na(coefs_ss_noint))] <- 0
  coefs_ss_int[which(is.na(coefs_ss_int))] <- 0
  coefs_ss_lasso[which(is.na(coefs_ss_lasso))] <- 0
  
  coefs_ss_noint_norm <- absnorm(coefs_ss_noint)
  
  outmat <- matrix(NA, ntest, 14)
  colnames(outmat) <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                        "Stack_noint", "Stack_noint_norm", "Stack_int",
                        "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  
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
    ss_lasso <- apply(allmod, 2, function(x){coefs_ss_lasso[1] + sum(coefs_ss_lasso[-1]*x)})
    ss_ridge <- apply(allmod, 2, function(x){coefs_ss_ridge[1] + sum(coefs_ss_ridge[-1]*x)})
    
    cury <- as.numeric(as.character(edat[[i]][,"y"]))
    
    outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2), mean((cury - unweighted)^2), mean((cury - sample_wtd)^2), mean((cury - cs_wtd)^2),
                                  mean((cury - stack_noint)^2), mean((cury - stack_noint_norm)^2), mean((cury - stack_int)^2),
                                  mean((cury - ss_noint)^2), mean((cury - ss_noint_norm)^2), 
                                  mean((cury - ss_int)^2), mean((cury - stack_lasso)^2), mean((cury - ss_lasso)^2), 
                                  mean((cury - stack_ridge)^2), mean((cury - ss_ridge)^2)))
  }
  outmat <- (outmat - outmat[,1])/outmat[,1]*100
  
  if (cluster_ind == 1){
    return(list(outmat = colMeans(outmat), index_max = index_max))
  }
  else{
    return(list(outmat = colMeans(outmat)))
  }
}


rep.clusters_fit <- function(reps, modfit, modpred, ndat, ntest, ncoef){
  #####
  logfile <- paste0("outputFile","nonsparsity",".txt")
  writeLines(c(""), file(logfile,'w'))
  
  num.threads <- as.integer(10)# round up
  threads <- makeCluster(num.threads, outfile=logfile)
  registerDoParallel(threads)
  
  getDoParWorkers()
  #######
  
  results <- foreach(i = 1:reps, .combine = 'rcomb', .multicombine = TRUE, .export = ls(globalenv())) %dopar% {
    library(locfit)
    library(fungible)
    library(mvtnorm)
    library(ClusterR)
    library(cluster)
    library(MLmetrics)
    library(randomForest)
    library(glmnet)
    library(nnet)
    library(nnls)
    library(plyr)
    library(BiocGenerics)
    library(clusterGeneration)
    sd <- sim_data(nstudies = ndat, ncoef = ncoef, ntest = ntest)$studies_list
    #for cluster: 
    cf <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 1, nsep = 0)
    errors_cluster <- cf$outmat
    indices <- cf$index_max+1
    errors_random<- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 2, nsep = indices)$outmat
    errors_multi <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 3, nsep = 0)$outmat
    return(list(errors_cluster = errors_cluster, errors_random = errors_random, errors_multi = errors_multi, indices = indices))
  }
  closeAllConnections()
  
  errors_cluster = na.omit(results[[1]])
  errors_random = na.omit(results[[2]])
  errors_multi = na.omit(results[[3]])
  indices <- as.vector(results[[4]])
  
  colnames(errors_multi) <- colnames(errors_cluster) <- colnames(errors_random) <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                                                                                     "Stack_noint", "Stack_noint_norm", "Stack_int",
                                                                                     "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  
  means_multi <- colMeans(errors_multi)
  means_cluster <- colMeans(errors_cluster)
  means_random <- colMeans(errors_random)
  sds_multi <- apply(errors_multi, 2, sd)
  sds_cluster <- apply(errors_cluster, 2, sd)
  sds_random <- apply(errors_random, 2, sd)
  
  return(list(means_multi = means_multi, means_cluster = means_cluster, means_random = means_random,
              sds_multi = sds_multi, sds_cluster = sds_cluster, sds_random = sds_random,
              errors_multi = errors_multi, errors_cluster = errors_cluster, errors_random = errors_random,
              indices = indices))   
}


vary_levels <- function(reps, var_list, modfit, modpred, ntest, out_str){
  ptm = proc.time()
  colnames_total <- c("Merged", "Unweighted", "Sample_Weighted", "CS_Weighted",
                      "Stack_noint", "Stack_noint_norm", "Stack_int",
                      "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  total_means_multi <- total_means_cluster <- total_means_random <- array(0, c(length(var_list), 14))
  total_sds_multi <- total_sds_cluster <- total_sds_random <- array(0, c(length(var_list), 14))
  
  colnames(total_means_multi) <- colnames(total_means_cluster) <- colnames(total_means_random) <- colnames_total 
  colnames(total_sds_multi) <- colnames(total_sds_cluster) <- colnames(total_sds_random) <- colnames_total
  
  indices_mat <- array(0, c(length(var_list), 2))
  colnames(indices_mat) <- c("mean", "sd")
  
  
  for (i in 1:length(var_list)){
    level <- var_list[i]
    print(level)
    
    level_rep <- rep.clusters_fit(reps, modfit, modpred, ndat = level, ntest, ncoef = 20)
    #means
    total_means_multi[i,] <- level_rep$means_multi
    total_means_cluster[i,] <- level_rep$means_cluster
    total_means_random[i,] <- level_rep$means_random
    #sds
    total_sds_multi[i,] <- level_rep$sds_multi
    total_sds_cluster[i,] <- level_rep$sds_cluster
    total_sds_random[i,] <- level_rep$sds_random
    indices_mat[i,] <- c(mean(level_rep$indices), sd(level_rep$indices))
    
    write.table(level_rep$errors_multi, paste0(out_str,"_errors_multi",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_cluster, paste0(out_str,"_errors_cluster",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_random, paste0(out_str,"_errors_random",level,".csv"), sep = ",", col.names = colnames_total)
  }
  
  
  colnames(total_means_multi) <- colnames(total_means_cluster) <- colnames(total_means_random) <- colnames_total 
  colnames(total_sds_multi) <- colnames(total_sds_cluster) <- colnames(total_sds_random) <- colnames_total
  
  write.table(total_means_multi, paste0(out_str,"_means_multi.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_cluster, paste0(out_str,"_means_cluster.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_random, paste0(out_str,"_means_random.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  write.table(total_sds_multi, paste0(out_str,"_sds_multi.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_cluster, paste0(out_str,"_sds_cluster.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_random, paste0(out_str,"_sds_random.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  write.table(indices_mat, paste0(out_str,"_numclusters.csv"), sep = ",", row.names = var_list, col.names = c("mean", "sd"))
  
  return(list(total_means_multi = total_means_multi, total_means_cluster = total_means_cluster, total_means_random = total_means_random,
              total_sds_multi = total_sds_multi, total_sds_cluster = total_sds_cluster, total_sds_random = total_sds_random,
              numclusters = indices_mat))
}


nclust_list <- c(7:10, 12, 14, 16, 18, 20, seq(30, 80, 10))
cs <- vary_levels(reps = 250, var_list = nclust_list, modfit = randomforestfit, modpred = randomforestpredict, ntest = 5, out_str = "cs")
