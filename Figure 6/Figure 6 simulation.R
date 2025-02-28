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
library(genefilter)
library(curatedOvarianData)
library(rpart)
library(Biobase)
library(dplyr)

#For home computer
#load("~/Desktop/edat_orig.Rda")

#For cluster
load("/home/mr356/edat_orig.Rda")


# Work with the intersection of rows
cn <- lapply(esets, rownames)
cn_int <- Reduce(intersect, cn)

edat <- vector("list", length(esets))

for(i in 1:length(esets)){
  edat[[i]] <- esets[[i]][cn_int,]
}

# Remove esets with missing gene expression data
ridx <- which(unlist(lapply(edat, function(x){sum(is.na(exprs(x)))})) > 0)
edat <- edat[-ridx] # length 15

# Convert eset list to set of matrices
for(i in 1:length(edat)){
  edat[[i]] <- t(exprs(edat[[i]]))
}

edat_orig <- edat

# Normalize the columns
for(i in 1:length(edat_orig)){
  edat_orig[[i]] <- apply(edat_orig[[i]], 2, scale)
}


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

randomforestfit <- function(data, ...){
  rf <- randomForest::randomForest(y ~ ., data = as.data.frame(data), ntree = 100, importance = TRUE, ...)
  rf
}

randomforestpredict <- function(data, newdata){
  as.vector(predict(data, newdata=as.data.frame(newdata)))
}

nnetfit <- function(data, ...){
  nnet::nnet(y ~ ., data = data, size = 10, MaxNWts = 2000, linout = T, trace = F)
}

nnetpred <- function(mod, newdata){
  as.vector(predict(mod, newdata = newdata))
}

ridgefit <- function(data, ...){
  mod <- glmnet::cv.glmnet(x = as.matrix(data[,-1]), y = data[,1], alpha = 0, ...)
}

ridgepred <- function(mod, newdata){
  as.vector(predict(mod, newx=as.matrix(newdata), s="lambda.1se"))
}


create_clusters <- function(studies_list, ntest, ncoef, k){
  merged <- do.call(rbind, studies_list[1:(length(studies_list) - ntest)])
  merged <- merged[sample(nrow(merged)), ]
  #cluster without using y
  #ss <- silhouette_score(k, merged[,-1])
  
  k2 <- kmeans(merged[,-1], centers = k, nstart = 25, iter.max = 25)
  
  if (k2$ifault==4) { 
    k2 = kmeans(merged[,-1], k2$centers, nstart = 25, iter.max = 25, algorithm="MacQueen")$cluster 
  }else{
    k2 <- k2$cluster
  }
  
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



#### Simulating outcomes and splitting into training and testing datasets 
sim_data <- function(edat_orig, ncoef, ntest){
  ndat = 15
  simtype = "nonl"
  #simtype = "normal"
  ninter = 10
  good = .25
  bad =  1
  val = .4
  icoefs <- c(4.4, 1.8, 2.5)
  setnc = TRUE 
  bin = FALSE
  
  edat <- edat_orig
  edat <- edat[sample(1:ndat)] 
  
  idx <- sample(1:ncol(edat[[1]]), ncoef)
  for(i in 1:ndat){
    edat[[i]] <- edat[[i]][,idx] 
  }
  
  if (setnc == TRUE){ #setting the number of coefficients used to generate outcome
    nchoose <- 10} #hardcoded: need to change this number
  else{#random number of coefficients used to generate the outcome
    #Generate Y - outcome for each patient given the genes selected 
    if(simtype == "nonl"){
      nchoose <- sample(3:ncoef, 1)
    } else {
      nchoose <- sample(2:ncoef, 1) #number of coefficients used to create the outcome is between 2 and nvar, sampled randomly
    }}
  
  coefs <- sample(c(runif(round(nchoose/2), -5, -0.5), runif(nchoose - round(nchoose/2), 0.5, 5))) 
  vars <- sample(1:ncol(edat[[1]]), nchoose) 
  
  #the mean percentage of the outcome (y) that is explained by the interaction terms
  for(i in 1:ndat){ #ndat is between 1 and 15 (number of datasets)
    
    curcoefs <- sapply(coefs, function(x){runif(1, x - (i<=5)*good - 
                                                  (i > 5 & i <= 10)*bad - (i > 10)*val , x + 
                                                  (i<=5)*good + (i > 5 & i <= 10)*bad + (i > 10)*val)}) #adds noise to the coefficients 
    
    
    if(simtype == "slash"){
      y <- (edat[[i]][,vars] %*% curcoefs) + 
        0.05*cbind(rnorm(nrow(edat[[i]]))/runif(nrow(edat[[i]]))) # Added "slash" noise
    } else if(simtype == "nonl"){
      if ((1 <= i) & (i <= ninter)){
        y <- (edat[[i]][,vars] %*% curcoefs) + icoefs[1]*edat[[i]][,vars[1]]*edat[[i]][,vars[2]] 
        - icoefs[2]*edat[[i]][,vars[1]]*edat[[i]][,vars[3]] + #(icoefs[3]*edat[[i]][,vars[3]])^2 + 
          10*sin(10*pi*edat[[i]][,vars[1]]) + cbind(rnorm(nrow(edat[[i]]))) # Added interaction terms
      }
      else{
        y <- (edat[[i]][,vars] %*% curcoefs) + cbind(rnorm(nrow(edat[[i]]))) # Added noise
      }
    } else {
      #y <- 10*sin(pi*edat[[i]][,2]* edat[[i]][,3])+20*(edat[[i]][,4] - .05)^2+10*edat[[i]][,5]+5*edat[[i]][,6] + cbind(rnorm(nrow(edat[[i]])))
      #y <- runif(1, 8, 12)*sin(pi*edat[[i]][,2]* edat[[i]][,3])+runif(1, 18, 22)*(edat[[i]][,4] - .05)^2+runif(1, 8, 12)*edat[[i]][,5]+runif(1, 4, 6)*edat[[i]][,6] + cbind(rnorm(nrow(edat[[i]])))
      
      y <- (edat[[i]][,vars] %*% curcoefs) + cbind(rnorm(nrow(edat[[i]]))) # Added noise
    }
    if (bin == TRUE){
      new_y <- ifelse(y >= quantile(y, .75), 1, 0)
      
      edat[[i]] <- cbind(new_y, edat[[i]])
      edat[[i]] <- as.data.frame(edat[[i]])
      colnames(edat[[i]]) <- c("y", paste0("V", 1:ncoef))
    }
    else{
      edat[[i]] <- cbind(y, edat[[i]])
      edat[[i]] <- as.data.frame(edat[[i]])
      colnames(edat[[i]]) <- c("y", paste0("V", 1:ncoef))
    }
  }
  if (bin == TRUE){
    for(i in 1:ndat){
      edat[[i]]$y <- factor(edat[[i]]$y)
    }
  }
  return(list(edat = edat, idx = idx, vars = vars, curcoefs = curcoefs))
}


#cluster_ind = 1 if running the algorithm on the k-means clusters, 2 if randomly generating substudies, 3 if just running on the original set of simulated studies
#rf_ind = 1 if the SSL is random forest, 2 if the SSL is ridge 
clusters_fit <- function(modfit, modpred, ndat, ncoef, ntest, studies_list, cluster_ind, rf_ind, nsep){
  #edat <- sim_data(ndat, ncoef)$studies_list
  if (cluster_ind == 1){ #K-means clustering
    #edat <- create_clusters(studies_list, ndat - ntest, ntest)$clusters_list
    cc <- create_clusters(studies_list, ntest, ncoef, nsep)
    edat <- cc$clusters_list
    #index_max <- cc$index_max
  }
  else if (cluster_ind == 2){ #random 'pseudo-studies'
    edat <- create_random(studies_list, ntest, ncoef, nsep)
  }
  else if (cluster_ind == 3){ #multistudy learning KNOWING the true clusters
    #edat <- studies_list
    edat <- studies_list
  }
  
  ntrain = length(edat) - ntest
  
  mods <- vector("list", ntrain)
  mses <- matrix(NA, ntrain, ntrain)
  
  allpreds <- vector("list", ntrain)
  
  #learn the merged predictor
  matstack <- do.call(rbind, edat[1:ntrain])
  #matstack <- matstack[sample(nrow(matstack)), ]
  
  #mod0 <- modfit(matstack[sample(nrow(matstack)), ])
  if (rf_ind == 1){
    mod0 <- randomForest::randomForest(y ~ ., data = as.data.frame(matstack[sample(nrow(matstack)), ]), ntree = 500, importance = TRUE)
  }
  else if (rf_ind == 2){
    mod0 <- modfit(matstack[sample(nrow(matstack)), ])
  }
  
  
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
  #outmat <- (outmat - outmat[,1])/outmat[,1]*100
  
  if (cluster_ind == 1){
    return(list(outmat = colMeans(outmat)))
  }
  else{
    return(list(outmat = colMeans(outmat)))
  }
  #return(outmat)
}

#test
#clusters_fit(modfit = ridgefit, modpred = ridgepred, ndat = 15, ncoef = 20, ntest = 5, studies_list = sd, cluster_ind = 1, rf_ind = 2, nsep = 5)

rep.clusters_fit <- function(reps, edat_orig, modfit, modpred, ndat, ntest, ncoef, k){
  #####
  logfile <- paste0("outputFile","nonsparsity",".txt")
  writeLines(c(""), file(logfile,'w'))
  
  #num.threads <- as.integer(reps/10)# round up
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
    sd <- sim_data(edat_orig, ncoef = ncoef, ntest = ntest)$edat
    #for cluster: 
    cf <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 1, rf_ind = 1, nsep = k)
    errors_cluster <- cf$outmat
    #indices <- c(indices, (cf$index_max+1))
    #indices <- cf$index_max+1
    errors_random<- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 2, rf_ind = 1, nsep = k)$outmat
    errors_multi <- clusters_fit(modfit, modpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 3, rf_ind = 1, nsep = k)$outmat
    errors_nnet <- clusters_fit(modfit = nnetfit, modpred = nnetpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 3, rf_ind = 2, nsep = k)$outmat
    #errors_ridge <- clusters_fit(modfit = ridgefit, modpred = ridgepred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 3, rf_ind = 2, nsep = k)$outmat
    errors_ridge <- clusters_fit(modfit = nnetfit, modpred = nnetpred, ndat, ncoef, ntest, studies_list = sd, cluster_ind = 1, rf_ind = 2, nsep = k)$outmat
    
    return(list(errors_cluster = errors_cluster, errors_random = errors_random, errors_multi = errors_multi, 
                errors_ridge = errors_ridge, errors_nnet = errors_nnet))
  }
  closeAllConnections()
  
  errors_cluster = (results[[1]])[,-3]
  errors_random = (results[[2]])[,-3]
  errors_multi = (results[[3]])[,-3]
  errors_ridge = (results[[4]])[,-3]
  errors_nnet = (results[[5]])[,-3]
  
  colnames(errors_multi) <- colnames(errors_cluster) <- colnames(errors_random) <- c("Merged", "Unweighted", "CS_Weighted",
                                                                                     "Stack_noint", "Stack_noint_norm", "Stack_int",
                                                                                     "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  
  means_multi <- colMeans(errors_multi)
  means_cluster <- colMeans(errors_cluster)
  means_random <- colMeans(errors_random)
  means_ridge <- colMeans(errors_ridge)
  means_nnet <- colMeans(errors_nnet)
  
  sds_multi <- apply(errors_multi, 2, sd)
  sds_cluster <- apply(errors_cluster, 2, sd)
  sds_random <- apply(errors_random, 2, sd)
  sds_ridge <- apply(errors_ridge, 2, sd)
  sds_nnet <- apply(errors_nnet, 2, sd)
  
  return(list(means_multi = means_multi, means_cluster = means_cluster, means_random = means_random, 
              means_ridge = means_ridge, means_nnet = means_nnet,
              sds_multi = sds_multi, sds_cluster = sds_cluster, sds_random = sds_random,
              sds_ridge = sds_ridge, sds_nnet = sds_nnet, 
              errors_multi = errors_multi, errors_cluster = errors_cluster, errors_random = errors_random,
              errors_ridge = errors_ridge, errors_nnet = errors_nnet))   
}


#t1 <- rep.clusters_fit(10, edat_orig, modfit = randomforestfit, modpred = randomforestpredict, ndat = 10, ntest = 5, ncoef = 20, k = 35)
#print(t1)

vary_levels <- function(reps, edat_orig, var_list, modfit, modpred, ndat, ntest, out_str){
  ptm = proc.time()
  colnames_total <- c("Merged", "Unweighted", "CS_Weighted",
                      "Stack_noint", "Stack_noint_norm", "Stack_int",
                      "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  total_means_multi <- total_means_cluster <- total_means_random <- total_means_ridge <- total_means_nnet <- array(0, c(length(var_list), 13))
  total_sds_multi <- total_sds_cluster <- total_sds_random <- total_sds_ridge <- total_sds_nnet <- array(0, c(length(var_list), 13))
  
  colnames(total_means_multi) <- colnames(total_means_cluster) <- colnames(total_means_random) <- 
    colnames(total_means_ridge) <- colnames(total_means_nnet) <- colnames_total 
  colnames(total_sds_multi) <- colnames(total_sds_cluster) <- colnames(total_sds_random) <- 
    colnames(total_sds_ridge) <- colnames(total_sds_nnet) <- colnames_total
  
  #indices_mat <- array(0, c(length(var_list), 2))
  #colnames(indices_mat) <- c("mean", "sd")
  
  for (i in 1:length(var_list)){
    level <- var_list[i]
    print(level)
    
    level_rep <- rep.clusters_fit(reps, edat_orig, modfit, modpred, ndat, ntest, ncoef = 20, k = level)
    #means
    total_means_multi[i,] <- level_rep$means_multi
    total_means_cluster[i,] <- level_rep$means_cluster
    total_means_random[i,] <- level_rep$means_random
    total_means_ridge[i, ] <- level_rep$means_ridge
    total_means_nnet[i, ] <- level_rep$means_nnet
    #sds
    total_sds_multi[i,] <- level_rep$sds_multi
    total_sds_cluster[i,] <- level_rep$sds_cluster
    total_sds_random[i,] <- level_rep$sds_random
    total_sds_ridge[i, ] <- level_rep$sds_ridge
    total_sds_nnet[i, ] <- level_rep$sds_nnet
    #indices
    #indices_mat[i,] <- c(mean(level_rep$indices), sd(level_rep$indices))
    
    write.table(level_rep$errors_multi, paste0(out_str,"_errors_multi",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_cluster, paste0(out_str,"_errors_cluster",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_random, paste0(out_str,"_errors_random",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_ridge, paste0(out_str,"_errors_ridge",level,".csv"), sep = ",", col.names = colnames_total)
    write.table(level_rep$errors_nnet, paste0(out_str,"_errors_nnet",level,".csv"), sep = ",", col.names = colnames_total)
  }
  write.table(total_means_multi, paste0(out_str,"_means_multi.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_cluster, paste0(out_str,"_means_cluster.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_random, paste0(out_str,"_means_random.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_ridge, paste0(out_str,"_means_ridge.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_means_nnet, paste0(out_str,"_means_nnet.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  write.table(total_sds_multi, paste0(out_str,"_sds_multi.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_cluster, paste0(out_str,"_sds_cluster.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_random, paste0(out_str,"_sds_random.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_ridge, paste0(out_str,"_sds_ridge.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  write.table(total_sds_nnet, paste0(out_str,"_sds_nnet.csv"), sep = ",", row.names = var_list, col.names = colnames_total)
  
  #write.table(indices_mat, paste0(out_str,"_numclusters.csv"), sep = ",", row.names = var_list, col.names = c("mean", "sd"))
  
  return(list(total_means_multi = total_means_multi, total_means_cluster = total_means_cluster, 
              total_means_random = total_means_random, total_means_ridge = total_means_ridge, total_means_nnet = total_means_nnet, 
              total_sds_multi = total_sds_multi, total_sds_cluster = total_sds_cluster, 
              total_sds_random = total_sds_random, total_sds_ridge = total_sds_ridge, total_sds_nnet = total_sds_nnet))
}




#Run the following:
k_list <- c(2, 5, 8, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100)
cs <- vary_levels(reps = 250, edat_orig, var_list = k_list, modfit = randomforestfit, modpred = randomforestpredict, ndat = 15, ntest = 5, out_str = "cs")

