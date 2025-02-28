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


# sim_data <- function(nstudies, ncoef, ntest){
#   studies_list  <- vector("list", nstudies)
#   nchoose <- 10
#   #nchoose <- ncoef
#   #general predictor-outcome rule:
#   coefs <- sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5)))
#   vars <- sample(1:ncoef, nchoose)
#   icoefs <- c(4, 1.8)
#   #If linear outcome model:
#   #icoefs <- c(0, 0)
# 
#   #norm.betas <- norm_vec(c(coefs,icoefs))
#   #coefs <- coefs/norm.betas
#   #cicoefs <- icoefs/norm.betas
# 
#   cormat <- matrix(.80, ncoef, ncoef)
#   diag(cormat) <- rep(1, ncoef)
#   in.cor.list <- replicate(nstudies, cormat, simplify=FALSE)
# 
# 
#   m <- monte(seed = sample(1:1000, 1), nvar = ncoef, nclus = nstudies,
#              clus.size = floor(runif(nstudies, 400, 400)), eta2 = runif(ncoef, .5, 1),
#              cor.list = NULL, random.cor = FALSE, skew.list = NULL,
#              kurt.list = NULL, secor = NULL, compactness = NULL,
#              sortMeans = FALSE)
# 
#   for (i in 1:nstudies){
#     curcoefs <- sapply(coefs, function(x){runif(1, x - .5, x + .5)})
# 
#     data_i <- scale(as.data.frame(m$data[m$data[,1] == i, ][,-1]))
# 
#     #To binarize:
# 
#     data_i <- sapply(1:ncol(data_i), function(x) {ifelse(data_i[,x] <= quantile(data_i[,x], runif(1, .2, .8)), 0, 1)})
# 
# 
# 
#     #generate outcome
#     #scaled here, but the original variables are left unscaled for the clustering step
#     y <- as.matrix((data_i[,vars]) %*% curcoefs) +
#       #icoefs[1]*(data_i[,vars[1]])^2 + icoefs[2]*(data_i[,vars[2]])^2 +
#      # + icoefs[1]*data_i[,vars[1]]*data_i[,vars[2]]
#     #- icoefs[2]*data_i[,vars[1]]*data_i[,vars[3]] +
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

sim_data <- function(nstudies, ncoef, ntest){
  studies_list  <- vector("list", nstudies)
  nchoose <- 10
  #nchoose <- ncoef
  #general predictor-outcome rule:
  coefs <- sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5)))
  vars <- sample(1:ncoef, nchoose)
  icoefs <- c(4, 1.8)
  n.noise <- 5
  #If linear outcome model:
  #icoefs <- c(0, 0)
  
  #norm.betas <- norm_vec(c(coefs,icoefs))
  #coefs <- coefs/norm.betas
  #cicoefs <- icoefs/norm.betas
  
  m <- genRandomClust(numClust = nstudies - ntest,
                      sepVal=0.8,
                      numNonNoisy=ncoef - n.noise,
                      numNoisy=n.noise,
                      numOutlier=100,
                      numReplicate=1,
                      fileName="test",
                      clustszind=2,
                      clustSizeEq=200,
                      rangeN=c(400,410),
                      clustSizes=NULL,
                      covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),
                      rangeVar=c(1, 10))
  
  m2 <- genRandomClust(numClust = 2,
                       sepVal=0.8,
                       numNonNoisy=ncoef - n.noise,
                       numNoisy=n.noise,
                       numOutlier=50,
                       numReplicate=ntest,
                       fileName="test",
                       clustszind=2,
                       clustSizeEq=50,
                       rangeN=c(500,600),
                       clustSizes=NULL,
                       covMethod=c("eigen", "onion", "c-vine", "unifcorrmat"),
                       rangeVar=c(1, 10))
  
  outliers <- scale(as.data.frame(m$datList)[which(m$memList[[1]] == 0), ])
  outliers_list <- vector("list", nstudies)
  
  d <- 1:nrow(outliers)
  x <- seq_along(d)
  d1 <- split(d, ceiling(x/(nrow(outliers)/nstudies)))
  for (i in 1:length(d1)){
    outliers_list[[i]] <- as.data.frame(outliers[d1[[i]],])
  }
  
  for (i in 1:nstudies){
    curcoefs <- sapply(coefs, function(x){runif(1, x - .5, x + .5)})
    
    if (i <= (nstudies - ntest)){
      data_i <- as.data.frame(m$datList)[which(m$memList[[1]] == i), ]
      data_i  <- scale(rbind(data_i, outliers_list[[i]]))
    }
    else{
      data_i <-scale(as.data.frame(m2$datList[[i - (nstudies - ntest)]]))
    }
    #data_i <- scale(as.data.frame(m$X[m$id == i, ]))
    
    #generate outcome
    #scaled here, but the original variables are left unscaled for the clustering step
    y <- as.matrix((data_i[,vars]) %*% curcoefs) +
      #icoefs[1]*(data_i[,vars[1]])^2 + icoefs[2]*(data_i[,vars[2]])^2 +
      # + icoefs[1]*data_i[,vars[1]]*data_i[,vars[2]]
      #- icoefs[2]*data_i[,vars[1]]*data_i[,vars[3]] +
      cbind(rnorm(nrow(data_i))) # Added noise
    
    studies_list[[i]] <- as.data.frame(cbind(y, data_i))
    colnames(studies_list[[i]]) <- c("y", paste0("V", 1:ncoef))
  }
  return(list(studies_list = studies_list))
}

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


create_clusters <- function(studies_list, ntest, ncoef, k){
  merged <- do.call(rbind, studies_list[1:(length(studies_list) - ntest)])
  merged <- merged[sample(nrow(merged)), ]
  #cluster without using y
  #ss <- silhouette_score(k, merged[,-1])
  
  k2 <- kmeans(merged[,-1], centers = k, nstart = 25)$cluster
  #k2 <- ss$km_ind
  #index_max <- ss$index_max
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
  return(list(clusters_list = clusters_list))
}


create_random <- function(studies_list, ntest, ncoef, nsep){
  ntr <- length(studies_list) - ntest
  edat <-  vector("list",(nsep + ntest))
  merged <- do.call(rbind, studies_list[1:ntr])
  merged <- merged[sample(nrow(merged)), ]
  
  
  d <- (1:nrow(merged))[sample(1:nrow(merged))]
  x <- seq_along(d)
  d1 <- split(d, ceiling(x/(nrow(merged)/nsep)))
  for (i in 1:length(d1)){
    edat[[i]] <- as.data.frame(merged[d1[[i]],])
  }
  for (t in 1:ntest){
    edat[[(nsep + t)]] <- as.data.frame(studies_list[[(length(studies_list) - ntest + t)]])
  }
  return(edat)
}


#cluster_ind = 1 if running the algorithm on the k-means clusters, 2 if randomly generating substudies, 3 if just running on the original set of simulated studies
clusters_fit <- function(modfit, modpred, ndat, ncoef, ntest, studies_list, cluster_ind, nsep){
  #edat <- sim_data(ndat, ncoef)$studies_list
  if (cluster_ind == 1){ #K-means clustering
    #edat <- create_clusters(studies_list, ndat - ntest, ntest)$clusters_list
    cc <- create_clusters(studies_list, ntest, ncoef, nsep)
    edat <- cc$clusters_list
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
  coefs_ss_lasso <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 1, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_ridge <- mapply(function(x,y){as.vector(coef(glmnet::cv.glmnet(x = y,
                                                                          y = as.numeric(as.character(x[,"y"])), alpha = 0, lower.limits = 0, intercept = T)))},edat[1:ntrain], allpreds, SIMPLIFY = FALSE)
  coefs_ss_lasso <- colMeans(do.call(rbind, coefs_ss_lasso), na.rm = T)
  coefs_ss_ridge <- colMeans(do.call(rbind, coefs_ss_ridge), na.rm = T)
  coefs_ss_noint[which(is.na(coefs_ss_noint))] <- 0
  coefs_ss_int[which(is.na(coefs_ss_int))] <- 0
  coefs_ss_lasso[which(is.na(coefs_ss_lasso))] <- 0
  coefs_ss_noint_norm <- absnorm(coefs_ss_noint)
  
  return(list(coefs = coefs_stack_ridge))
}

rep.clusters_fit <- function(reps, modfit, modpred, ndat, ntest, ncoef, k){
  #####
  logfile <- paste0("outputFile","nonsparsity",".txt")
  writeLines(c(""), file(logfile,'w'))
  
  num.threads <- 10
  threads <- makeCluster(num.threads, outfile=logfile,setup_timeout = 0.5)
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
    errors_cluster <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 1, nsep = k)
    errors_random<- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 2, nsep = k)
    errors_multi <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 3, nsep = k)
    
    
    stack_cluster <- errors_cluster$coefs
    stack_random <- errors_random$coefs
    stack_multi <- errors_multi$coefs
    return(list(stack_cluster, stack_random, stack_multi))
  }
  closeAllConnections()
  
  stack_cluster <- results[[1]]
  stack_random <- results[[2]]
  stack_multi <- results[[3]]
  
  out_str <- "cs"
  
  write.table(stack_cluster, paste0(out_str,"_stack_cluster_",k,".csv"), sep = ",")
  write.table(stack_random, paste0(out_str,"_stack_random_",k,".csv"), sep = ",")
  write.table(stack_multi, paste0(out_str,"_stack_multi_",k,".csv"), sep = ",")
  
  return(list(stack_cluster, stack_random, stack_multi))   
}

#K = 20
rcf1 <- rep.clusters_fit(reps = 100, modfit = randomforestfit, modpred = randomforestpredict, ndat = 10, ntest = 5, ncoef = 20, k = 20)
#K = 80
rcf2 <- rep.clusters_fit(reps = 100, modfit = randomforestfit, modpred = randomforestpredict, ndat = 10, ntest = 5, ncoef = 20, k = 80)