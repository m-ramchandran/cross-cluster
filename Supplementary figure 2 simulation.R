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

randomforestpredict <- function(data, newdata, treeweight){
  #newdata = as.data.frame(newdata)
  if (treeweight == TRUE){
    pred_obj <- predict(data, newdata=newdata, predict.all = TRUE)
    pred_obj$individual
  }
  else{
    as.vector(predict(data, newdata=newdata))
  }
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

  cormat <- matrix(.80, ncoef, ncoef)
  diag(cormat) <- rep(1, ncoef)
  in.cor.list <- replicate(nstudies, cormat, simplify=FALSE)
  
  
  m <- monte(seed = sample(1:1000, 1), nvar = ncoef, nclus = nstudies,
             clus.size = floor(runif(nstudies, 500, 510)), eta2 = runif(ncoef, .5, 1),
             cor.list = NULL, random.cor = FALSE, skew.list = NULL,
             kurt.list = NULL, secor = NULL, compactness = NULL,
             sortMeans = FALSE)
  
  for (i in 1:nstudies){
    curcoefs <- sapply(coefs, function(x){runif(1, x - .5, x + .5)})
    data_i <- scale(as.data.frame(m$data[m$data[,1] == i, ][,-1]))
    
    #To binarize:
    #data_i <- sapply(1:ncol(data_i), function(x) {ifelse(data_i[,x] <= quantile(data_i[,x], runif(1, .2, .8)), 0, 1)})
    
    #generate outcome
    #scaled here, but the original variables are left unscaled for the clustering step
    y <- as.matrix((data_i[,vars]) %*% curcoefs) +
      icoefs[1]*(data_i[,vars[1]])^2 + icoefs[2]*(data_i[,vars[2]])^2 +
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
  k2 <- kmeans(merged[,-1], centers = k, nstart = 25)$cluster
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
clusters_fit <- function(modfit, modpred, ndat, ncoef, ntest, studies_list, cluster_ind, nsep, treeweight){
  if (cluster_ind == 1){ #K-means clustering
    cc <- create_clusters(studies_list, ntest, ncoef, nsep)
    edat <- cc$clusters_list
  }
  else if (cluster_ind == 2){ #random 'pseudo-studies'
    edat <- create_random(studies_list, ntest, ncoef, nsep)
  }
  else if (cluster_ind == 3){ #multistudy learning KNOWING the true clusters
    edat <- studies_list
  }
  
  numtree <- 100
  ntrain = length(edat) - ntest
  mods <- vector("list", ntrain)
  
  if (treeweight == TRUE){
    mses <- matrix(NA, ntrain*numtree, ntrain)
  }
  else{
    mses <- matrix(NA, ntrain, ntrain)
  }
  
  allpreds <- vector("list", ntrain)
  
  #learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain])
  mod0 <- randomForest::randomForest(y ~ ., data = as.data.frame(matstack[sample(nrow(matstack)), ]), ntree = 500, importance = TRUE)
  
  
  if (treeweight == TRUE){
    for(i in 1: ntrain){
      mods[[i]] <- modfit(edat[[i]])
      preds <- lapply(edat[1:ntrain], function(x){
        modpred(mods[[i]], newdata = x[, -1], treeweight) 
      })
      new_pred <- vector("list", numtree)
      for (j in 1:numtree){
        new_pred[[j]] <- vector("list", ntrain)
      }
      for (a in 1:ntrain){
        preds_k <- preds[[a]]
        for (h in 1:numtree){
          new_pred[[h]][[a]] <- as.vector(preds_k[,h])
        }
      }
      
      for (j in 1:numtree){
        if (i == 1){p <- j}
        else{p <- (i - 1)*numtree + j}
        mses[p,] <- unlist(lapply(1:ntrain, function(i){
          x = edat[[i]]
          mean((new_pred[[j]][[i]] - x[,"y"])^2)
        }))
        curpreds <- lapply(new_pred[[j]], as.numeric)
        allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)
      }}
    for (i in 1:ntrain){
      h <- (i-1)*numtree + 1 
      mses[,i][h:(h+numtree-1)] <- NA 
    }}
  else{
    for(i in 1:ntrain){
      mods[[i]] <- modfit(edat[[i]])
      preds <- lapply(edat[1:ntrain], function(x){
        modpred(mods[[i]], newdata = x[, -1], treeweight) 
      })
      mses[i,] <- unlist(lapply(edat[1:ntrain], function(x){#cross validation within the training set
        newdata = x[, -1]
        preds <- modpred(mods[[i]], newdata = newdata, treeweight) 
        mean((preds - x[,"y"])^2)}
      ))
      
      curpreds <- lapply(preds, as.numeric)
      allpreds <- mapply(cbind, allpreds, curpreds, SIMPLIFY = FALSE)}
    diag(mses) <- NA
  }
  
  # CS Weights
  tt <- apply(mses, 1, mean, na.rm = T) #removing the diagonal elements, takes the mean of each row 
  weights <- absnorm(sqrt(tt), max.norm = TRUE)
  nk <- unlist(lapply(edat, nrow)) #vector of number of rows in each dataset
  nwts <- absnorm(nk[1:ntrain])
  
  # Regression: stacked (intercept and no intercept)
  predstack <- do.call(rbind, allpreds)
  
  #coefs_stack_noint <- nnls::nnls(predstack, as.numeric(as.character(matstack$y)))$x
  coefs_stack_int <- nnls::nnls(cbind(rep(1,nrow(predstack)),predstack), as.numeric(as.character(matstack$y)))$x
  coefs_stack_lasso <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 1, lower.limits = 0, intercept = T)))
  coefs_stack_ridge <- as.vector(coef(glmnet::cv.glmnet(x = predstack, y = as.numeric(as.character(matstack$y)), alpha = 0, lower.limits = 0, intercept = T)))
  #Just a safeguard against full collinearity, although I think we are OK with nnls now
  coefs_stack_int[which(is.na(coefs_stack_int))] <- 0
  coefs_stack_lasso[which(is.na(coefs_stack_lasso))] <- 0
  coefs_stack_ridge[which(is.na(coefs_stack_ridge))] <- 0

  outmat <- matrix(NA, ntest, 4)
  colnames(outmat) <- c("Merged", "Stack_int","Stack_lasso","Stack_ridge")

  for(i in (ntrain + 1):(length(edat))){
    merged <- modpred(mod0, newdata = edat[[i]][,-1], treeweight)
    if (treeweight == TRUE){
      merged <- apply(merged, 2, as.numeric)
      allmod <- t(do.call(cbind,lapply(mods, modpred, newdata = edat[[i]][,-1], treeweight)))
      allmod <- apply(allmod, 2, as.numeric)
    }
    else{
      merged <- as.vector(sapply(merged, as.numeric))
      allmod <- do.call(rbind,lapply(mods, modpred, newdata = edat[[i]][,-1], treeweight))
      allmod <- apply(allmod, 2, as.numeric)
    }
    stack_int <- apply(allmod, 2, function(x){coefs_stack_int[1] + sum(coefs_stack_int[-1]*x)})
    stack_lasso <- apply(allmod, 2, function(x){coefs_stack_lasso[1] + sum(coefs_stack_lasso[-1]*x)})
    stack_ridge <- apply(allmod, 2, function(x){coefs_stack_ridge[1] + sum(coefs_stack_ridge[-1]*x)})

    cury <- as.numeric(as.character(edat[[i]][,"y"]))
    outmat[i - ntrain,] <- sqrt(c(mean((cury - merged)^2), mean((cury - stack_int)^2),
                                  mean((cury - stack_lasso)^2), 
                                  mean((cury - stack_ridge)^2)))
  }
  outmat <- (outmat - outmat[,1])/outmat[,1]*100
  
  if (cluster_ind == 1){
    return(list(outmat = colMeans(outmat)))
  }
  else{
    return(list(outmat = colMeans(outmat)))
  }
}

rep.clusters_fit <- function(reps, modfit, modpred, ndat, ntest, ncoef, k){
  #####
  logfile <- paste0("outputFile","nonsparsity",".txt")
  writeLines(c(""), file(logfile,'w'))
  
  #num.threads <- round(reps/10)# round up
  num.threads <- as.integer(10)
  threads <- makeCluster(num.threads, outfile=logfile, setup_timeout = 0.5)
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
    #Weighting forests
    errors_cluster_forest <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 1, nsep = k, treeweight = FALSE)$outmat
    errors_random_forest<- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 2, nsep = k, treeweight = FALSE)$outmat
    errors_multi_forest <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 3, nsep = k, treeweight = FALSE)$outmat
    
    
    #Weighting trees 
    errors_cluster_tree <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 1, nsep = k, treeweight = TRUE)$outmat
    errors_random_tree<- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 2, nsep = k, treeweight = TRUE)$outmat
    errors_multi_tree <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 3, nsep = k, treeweight = TRUE)$outmat
    return(list(errors_cluster_forest, errors_random_forest, errors_multi_forest, errors_cluster_tree, errors_random_tree, errors_multi_tree))
  }
  closeAllConnections()

  errors_cluster_forest = results[[1]]
  errors_random_forest = results[[2]]
  errors_multi_forest = results[[3]]
  
  errors_cluster_tree = results[[4]]
  errors_random_tree = results[[5]]
  errors_multi_tree = results[[6]]

  colnames_total <- c("Merged", "Stack_int","Stack_lasso","Stack_ridge")
  
  
  colnames(errors_multi_forest) <- colnames(errors_cluster_forest) <- colnames(errors_random_forest) <-
    colnames(errors_multi_tree) <- colnames(errors_cluster_tree) <- colnames(errors_random_tree) <- colnames_total
  
  #Forest
  means_multi_forest <- colMeans(errors_multi_forest)
  means_cluster_forest <- colMeans(errors_cluster_forest)
  means_random_forest <- colMeans(errors_random_forest)
  sds_multi_forest <- apply(errors_multi_forest, 2, sd)
  sds_cluster_forest <- apply(errors_cluster_forest, 2, sd)
  sds_random_forest <- apply(errors_random_forest, 2, sd)
  
  #Tree
  means_multi_tree <- colMeans(errors_multi_tree)
  means_cluster_tree <- colMeans(errors_cluster_tree)
  means_random_tree <- colMeans(errors_random_tree)
  sds_multi_tree <- apply(errors_multi_tree, 2, sd)
  sds_cluster_tree <- apply(errors_cluster_tree, 2, sd)
  sds_random_tree <- apply(errors_random_tree, 2, sd)
  
  return(list(means_multi_forest = means_multi_forest, means_cluster_forest = means_cluster_forest, means_random_forest = means_random_forest,
              sds_multi_forest = sds_multi_forest, sds_cluster_forest = sds_cluster_forest, sds_random_forest = sds_random_forest,
              means_multi_tree = means_multi_tree, means_cluster_tree = means_cluster_tree, means_random_tree = means_random_tree,
              sds_multi_tree = sds_multi_tree, sds_cluster_tree = sds_cluster_tree, sds_random_tree = sds_random_tree,
              errors_multi_forest = errors_multi_forest, errors_cluster_forest = errors_cluster_forest, errors_random_forest = errors_random_forest,
              errors_multi_tree = errors_multi_tree, errors_cluster_tree = errors_cluster_tree, errors_random_tree = errors_random_tree))   
}

vary_levels <- function(reps, var_list, modfit, modpred, ndat, ntest, out_str){
  ptm = proc.time()
  
  colnames_total <- c("Merged", "Stack_int","Stack_lasso","Stack_ridge")

  total_means_multi_forest <- total_means_cluster_forest <- total_means_random_forest <- 
    total_sds_multi_forest <- total_sds_cluster_forest <- total_sds_random_forest <- 
    total_means_multi_tree <- total_means_cluster_tree <- total_means_random_tree <- 
    total_sds_multi_tree <- total_sds_cluster_tree <- total_sds_random_tree <- array(0, c(length(var_list), 4))
  
  colnames(total_means_multi_forest) <- colnames(total_means_cluster_forest) <- colnames(total_means_random_forest) <- colnames_total 
  colnames(total_sds_multi_forest) <- colnames(total_sds_cluster_forest) <- colnames(total_sds_random_forest) <- colnames_total
  colnames(total_means_multi_tree) <- colnames(total_means_cluster_tree) <- colnames(total_means_random_tree) <- colnames_total 
  colnames(total_sds_multi_tree) <- colnames(total_sds_cluster_tree) <- colnames(total_sds_random_tree) <- colnames_total
  
  for (i in 1:length(var_list)){
    level <- var_list[i]
    print(level)
    
    level_rep <- rep.clusters_fit(reps, modfit, modpred, ndat, ntest, ncoef = 20, k = level)
    #Forest means
    total_means_multi_forest[i,] <- level_rep$means_multi_forest
    total_means_cluster_forest[i,] <- level_rep$means_cluster_forest
    total_means_random_forest[i,] <- level_rep$means_random_forest
    #Forest sds
    total_sds_multi_forest[i,] <- level_rep$sds_multi_forest
    total_sds_cluster_forest[i,] <- level_rep$sds_cluster_forest
    total_sds_random_forest[i,] <- level_rep$sds_random_forest
    
    #Tree means
    total_means_multi_tree[i,] <- level_rep$means_multi_tree
    total_means_cluster_tree[i,] <- level_rep$means_cluster_tree
    total_means_random_tree[i,] <- level_rep$means_random_tree
    #Tree sds
    total_sds_multi_tree[i,] <- level_rep$sds_multi_tree
    total_sds_cluster_tree[i,] <- level_rep$sds_cluster_tree
    total_sds_random_tree[i,] <- level_rep$sds_random_tree
    
    #Forest errors
    write.table(level_rep$errors_multi_forest, paste0(out_str,"_errors_multi_forest",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_cluster_forest, paste0(out_str,"_errors_cluster_forest",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_random_forest, paste0(out_str,"_errors_random_forest",level,".csv"), sep = ",", col.names = colnames_total)
    #Tree errors
    write.table(level_rep$errors_multi_tree, paste0(out_str,"_errors_multi_tree",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_cluster_tree, paste0(out_str,"_errors_cluster_tree",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_random_tree, paste0(out_str,"_errors_random_tree",level,".csv"), sep = ",", col.names = colnames_total)
  }
  
  colnames(total_means_multi_forest) <- colnames(total_means_cluster_forest) <- colnames(total_means_random_forest) <- colnames_total 
  colnames(total_sds_multi_forest) <- colnames(total_sds_cluster_forest) <- colnames(total_sds_random_forest) <- colnames_total
  colnames(total_means_multi_tree) <- colnames(total_means_cluster_tree) <- colnames(total_means_random_tree) <- colnames_total 
  colnames(total_sds_multi_tree) <- colnames(total_sds_cluster_tree) <- colnames(total_sds_random_tree) <- colnames_total
  
  #Forest
  write.table(total_means_multi_forest, paste0(out_str,"_means_multi_forest.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_cluster_forest, paste0(out_str,"_means_cluster_forest.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_random_forest, paste0(out_str,"_means_random_forest.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  write.table(total_sds_multi_forest, paste0(out_str,"_sds_multi_forest.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_cluster_forest, paste0(out_str,"_sds_cluster_forest.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_random_forest, paste0(out_str,"_sds_random_forest.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  #Tree
  write.table(total_means_multi_tree, paste0(out_str,"_means_multi_tree.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_cluster_tree, paste0(out_str,"_means_cluster_tree.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_random_tree, paste0(out_str,"_means_random_tree.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  write.table(total_sds_multi_tree, paste0(out_str,"_sds_multi_tree.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_cluster_tree, paste0(out_str,"_sds_cluster_tree.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_random_tree, paste0(out_str,"_sds_random_tree.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  return(list(total_means_multi_forest = total_means_multi_forest, total_means_cluster_forest = total_means_cluster_forest, total_means_random_forest = total_means_random_forest,
              total_sds_multi_forest = total_sds_multi_forest, total_sds_cluster_forest = total_sds_cluster_forest, total_sds_random_forest = total_sds_random_forest,
              total_means_multi_tree = total_means_multi_tree, total_means_cluster_tree = total_means_cluster_tree, total_means_random_tree = total_means_random_tree,
              total_sds_multi_tree = total_sds_multi_tree, total_sds_cluster_tree = total_sds_cluster_tree, total_sds_random_tree = total_sds_random_tree))
}



k_list <- c(2, 5, 8, seq(10, 80, 10))
cs <- vary_levels(reps = 100, var_list = k_list, modfit = randomforestfit, modpred = randomforestpredict, ndat = 10, ntest = 5, out_str = "cs")
