##############################################################
#Varying number of clusters
##############################################################
multi.plot1 <- function(ncoef_list, filename, title, xlab, ylab,xlim, ylim, n, ncoef_ind){
  paste0("~/Desktop/", filename, "/cs_means_cluster.csv")
  
  means_cluster <- read.csv(paste0("~/Desktop/", filename, "/cs_means_cluster.csv"))[, -1]
  sds_cluster <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_cluster.csv"))[, -1]
  means_multi <- read.csv(paste0("~/Desktop/", filename, "/cs_means_multi.csv"))[, -1]
  sds_multi <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_multi.csv"))[, -1]
  means_random <- read.csv(paste0("~/Desktop/", filename, "/cs_means_random.csv"))[, -1]
  sds_random <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_random.csv"))[, -1]
  
  #Creating the plotting data frames
  means_mat <- as.matrix(cbind(means_cluster$Stack_ridge, means_random$Stack_ridge, means_multi$Stack_ridge))
  if (ncoef_ind == TRUE){
    ncoef_list <- ncoef_list - 5
  }
  data_means <- data.frame(cbind(ncoef_list, means_mat))
  df.melted <- melt(data_means, id = "ncoef_list")
  se.cluster <- data.frame(cbind(ncoef_list, means_cluster$Stack_ridge - 1.96*sds_cluster$Stack_ridge/(sqrt(n)),
                                 means_cluster$Stack_ridge + 1.96*sds_cluster$Stack_ridge/(sqrt(n))))
  se.multi <- data.frame(cbind(ncoef_list, means_multi$Stack_ridge - 1.96*sds_multi$Stack_ridge/(sqrt(n)),
                               means_multi$Stack_ridge + 1.96*sds_multi$Stack_ridge/(sqrt(n))))
  se.random <- data.frame(cbind(ncoef_list, means_random$Stack_ridge - 1.96*sds_random$Stack_ridge/(sqrt(n)),
                                means_random$Stack_ridge + 1.96*sds_random$Stack_ridge/(sqrt(n))))
  
  #Plot
  q1 <- ggplot(data = df.melted, aes(x = ncoef_list, y = value, color = variable)) +
    geom_point() + geom_line() + 
    geom_ribbon(data=se.cluster,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#2A9D8F",alpha=0.2, inherit.aes = FALSE) + 
    geom_ribbon(data=se.random,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#F4A261",alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data=se.multi,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#254653",alpha=0.2, inherit.aes = FALSE) +
    theme_classic() + scale_y_continuous("", limits = ylim, breaks  = scales::pretty_breaks(n = 5)) + 
    scale_x_continuous(xlab, limits = xlim, breaks  = scales::pretty_breaks(n = 5)) + labs(title = title) +
    theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)),
          plot.title = element_text(size = rel(1.3)), legend.text=element_text(size=rel(1.3)), 
          legend.title=element_text(size=rel(1.5)), legend.position = "none") +
    scale_color_manual(values=c("#2A9D8F","#F4A261", "#254653"), labels = c("Cluster", "Random", "Multi"), name="Method")
  return(q1)
}
nclust_list <- c(7:10, 12, 14, 16, 18, 20, seq(30, 80, 10))
numclust <- multi.plot1(ncoef_list = nclust_list, filename = "250reps/numclust", title = "Varying number of true clusters", 
                        xlab = "Number of clusters", ylab = "% change in RMSE from Merged", xlim = c(0, 80), ylim = c(-40, 0),
                        n = 250, ncoef_ind = TRUE)
show(numclust)

##############################################################
#Varying signal
##############################################################
norm_list <- c(0, .5, 1, 1.5, 2, 3, 5, 7.5, 10, 15, 20, 25, 30)
signal <- multi.plot1(ncoef_list = norm_list, filename = "250reps/signal", title = "Increasing strength of signal", 
                      xlab = "Norm of coefficient vector", ylab = "% change in RMSE from Merged", xlim = c(0, 30), ylim = c(-60, 5),n = 250, ncoef_ind = FALSE)
show(signal)


##############################################################
#Varying sample size 
##############################################################
multi.plot.files <- function(ncoef_list, filename, title, xlab, ylab, xlim, ylim, n){
  
  means_cluster <- means_multi <- means_random <- as.data.frame(array(0, c(length(ncoef_list), 13)))
  sds_cluster <- sds_multi <- sds_random <- as.data.frame(array(0, c(length(ncoef_list), 13)))
  colnames_total <- c("Unweighted", "Sample_Weighted","CS_Weighted",
                      "Stack_noint", "Stack_noint_norm", "Stack_int",
                      "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  colnames(means_cluster) <- colnames(means_multi) <- colnames(means_random) <- colnames_total
  colnames(sds_cluster) <- colnames(sds_multi) <- colnames(sds_random) <- colnames_total
  
  for (i in (1:length(ncoef_list))){
    n <- ncoef_list[i]
    means_cluster[i,] <- apply(na.omit(read.csv(paste0("~/Desktop/", filename, "/cs_errors_cluster", n,".csv"))), 2, mean)[-1]
    means_multi[i,] <- apply(na.omit(read.csv(paste0("~/Desktop/", filename, "/cs_errors_multi", n,".csv"))), 2, mean)[-1]
    means_random[i,] <- apply(na.omit(read.csv(paste0("~/Desktop/", filename, "/cs_errors_random", n,".csv"))), 2, mean)[-1]
    sds_cluster[i,] <- apply(na.omit(read.csv(paste0("~/Desktop/", filename, "/cs_errors_cluster", n,".csv"))), 2, sd)[-1]
    sds_multi[i,] <- apply(na.omit(read.csv(paste0("~/Desktop/", filename, "/cs_errors_multi", n,".csv"))), 2, sd)[-1]
    sds_random[i,] <- apply(na.omit(read.csv(paste0("~/Desktop/", filename, "/cs_errors_random", n,".csv"))), 2, sd)[-1]
  }
  
  #Creating the plotting data frames
  means_mat <- as.matrix(cbind(means_cluster$Stack_ridge, means_random$Stack_ridge, means_multi$Stack_ridge))
  data_means <- data.frame(cbind(ncoef_list, means_mat))
  df.melted <- melt(data_means, id = "ncoef_list")
  se.cluster <- data.frame(cbind(ncoef_list, means_cluster$Stack_ridge - 1.96*sds_cluster$Stack_ridge/(sqrt(100)),
                                 means_cluster$Stack_ridge + 1.96*sds_cluster$Stack_ridge/(sqrt(100))))
  se.multi <- data.frame(cbind(ncoef_list, means_multi$Stack_ridge - 1.96*sds_multi$Stack_ridge/(sqrt(100)),
                               means_multi$Stack_ridge + 1.96*sds_multi$Stack_ridge/(sqrt(100))))
  se.random <- data.frame(cbind(ncoef_list, means_random$Stack_ridge - 1.96*sds_random$Stack_ridge/(sqrt(100)),
                                means_random$Stack_ridge + 1.96*sds_random$Stack_ridge/(sqrt(100))))
  
  #Plot
  q1 <- ggplot(data = df.melted, aes(x = ncoef_list, y = value, color = variable)) +
    geom_point() + geom_line() + 
    geom_ribbon(data=se.cluster,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#2A9D8F",alpha=0.2, inherit.aes = FALSE) + 
    geom_ribbon(data=se.random,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#F4A261",alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data=se.multi,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#254653",alpha=0.2, inherit.aes = FALSE) +
    theme_classic() + scale_y_continuous("", limits = ylim, breaks  = scales::pretty_breaks(n = 5)) + 
    scale_x_continuous(xlab, limits = xlim, breaks  = scales::pretty_breaks(n = 5)) + labs(title = title) +
    theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)),
          plot.title = element_text(size = rel(1.3)), legend.text=element_text(size=rel(1.3)), 
          legend.title=element_text(size=rel(1.5)), legend.position = "none") +
    scale_color_manual(values=c("#2A9D8F","#F4A261", "#254653"), labels = c("Cluster", "Random", "Multi"), name="Method")
  #scale_color_manual(name = "Method", values = c("Cluster" = "#2A9D8F",  "Random" = "#F4A261", "Multi" =  "#254653"))
  return(q1)
}

nsize_list <- seq(100, 1000, 100)
sampsize <- multi.plot.files(ncoef_list = nsize_list, filename = "250reps/sampsize", title = "Increasing sample size per cluster", 
                             xlab = "Sample size of each cluster", ylab = "% change in RMSE from Merged", xlim = c(0, 1000), ylim = c(-40, 0),n = 250)
show(sampsize)

##############################################################
#Composite graph for the varying parameters
##############################################################
library(ggpubr)
vary_params <- ggarrange(signal, numclust, sampsize, 
                         labels = c("A", "B", "C"),
                         ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
show(vary_params)

annotate_figure(vary_params,left = text_grob("% change in average RMSE from Merged", rot = 90, size = 14),
                top = text_grob("Varying dataset characteristics", size = 20, face = "bold"))


