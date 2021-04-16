library(ggplot2)
library(reshape2)
library(ggpubr)

##############################################################
#Treeweighting plot 
##############################################################
multi.plot.tree <- function(ncoef_list, filename, title, xlab, ylab,xlim, ylim, n){
  paste0("~/Desktop/", filename, "/cs_means_cluster.csv")
  
  colnames_total <- c("Merged", "Stack_int", "Stack_lasso","Stack_ridge")
  
  #Reading in the tree files 
  means_cluster_tree <- read.csv(paste0("~/Desktop/", filename, "/cs_means_cluster_tree.csv"))
  sds_cluster_tree <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_cluster_tree.csv"))
  
  means_multi_tree <- read.csv(paste0("~/Desktop/", filename, "/cs_means_multi_tree.csv"))
  sds_multi_tree <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_multi_tree.csv"))
  
  means_multi_tree <- do.call("rbind", replicate(nrow(means_multi_tree), means_multi_tree[5,], simplify = FALSE))
  sds_multi_tree <- do.call("rbind", replicate(nrow(sds_multi_tree), sds_multi_tree[5,], simplify = FALSE))
  
  means_random_tree <- read.csv(paste0("~/Desktop/", filename, "/cs_means_random_tree.csv"))
  sds_random_tree <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_random_tree.csv"))
  
  colnames(means_cluster_tree) <- colnames(means_multi_tree) <- colnames(means_random_tree) <- colnames_total
  colnames(sds_cluster_tree) <- colnames(sds_multi_tree) <- colnames(sds_random_tree) <- colnames_total
  
  #Reading in the forest files 
  means_cluster_forest <- read.csv(paste0("~/Desktop/", filename, "/cs_means_cluster_forest.csv"))
  sds_cluster_forest <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_cluster_forest.csv"))
  
  means_multi_forest <- read.csv(paste0("~/Desktop/", filename, "/cs_means_multi_forest.csv"))
  sds_multi_forest <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_multi_forest.csv"))
  
  means_multi_forest <- do.call("rbind", replicate(nrow(means_multi_forest), means_multi_forest[5,], simplify = FALSE))
  sds_multi_forest <- do.call("rbind", replicate(nrow(sds_multi_forest), sds_multi_forest[5,], simplify = FALSE))
  
  means_random_forest <- read.csv(paste0("~/Desktop/", filename, "/cs_means_random_forest.csv"))
  sds_random_forest <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_random_forest.csv"))
  
  colnames(means_cluster_forest) <- colnames(means_multi_forest) <- colnames(means_random_forest) <- colnames_total
  colnames(sds_cluster_forest) <- colnames(sds_multi_forest) <- colnames(sds_random_forest) <- colnames_total
  
  
  #Creating the plotting data frames
  means_mat <- as.matrix(cbind(means_cluster_tree$Stack_ridge, means_random_tree$Stack_ridge, means_multi_tree$Stack_ridge,
                               means_cluster_forest$Stack_ridge, means_random_forest$Stack_ridge, means_multi_forest$Stack_ridge))
  data_means <- data.frame(cbind(ncoef_list, means_mat))
  df.melted <- melt(data_means, id = "ncoef_list")
  se.cluster_tree <- data.frame(cbind(ncoef_list, means_cluster_tree$Stack_ridge - 1.96*sds_cluster_tree$Stack_ridge/(sqrt(n)),
                                      means_cluster_tree$Stack_ridge + 1.96*sds_cluster_tree$Stack_ridge/(sqrt(n))))
  se.multi_tree <- data.frame(cbind(ncoef_list, means_multi_tree$Stack_ridge - 1.96*sds_multi_tree$Stack_ridge/(sqrt(n)),
                                    means_multi_tree$Stack_ridge + 1.96*sds_multi_tree$Stack_ridge/(sqrt(n))))
  se.random_tree <- data.frame(cbind(ncoef_list, means_random_tree$Stack_ridge - 1.96*sds_random_tree$Stack_ridge/(sqrt(n)),
                                     means_random_tree$Stack_ridge + 1.96*sds_random_tree$Stack_ridge/(sqrt(n))))
  
  se.cluster_forest <- data.frame(cbind(ncoef_list, means_cluster_forest$Stack_ridge - 1.96*sds_cluster_forest$Stack_ridge/(sqrt(n)),
                                        means_cluster_forest$Stack_ridge + 1.96*sds_cluster_forest$Stack_ridge/(sqrt(n))))
  se.multi_forest <- data.frame(cbind(ncoef_list, means_multi_forest$Stack_ridge - 1.96*sds_multi_forest$Stack_ridge/(sqrt(n)),
                                      means_multi_forest$Stack_ridge + 1.96*sds_multi_forest$Stack_ridge/(sqrt(n))))
  se.random_forest <- data.frame(cbind(ncoef_list, means_random_forest$Stack_ridge - 1.96*sds_random_forest$Stack_ridge/(sqrt(n)),
                                       means_random_forest$Stack_ridge + 1.96*sds_random_forest$Stack_ridge/(sqrt(n))))
  
  
  means_mat1 <- as.matrix(cbind(means_cluster_tree$Stack_ridge, means_random_tree$Stack_ridge,
                                means_cluster_forest$Stack_ridge, means_random_forest$Stack_ridge))
  data_means1 <- data.frame(cbind(ncoef_list, means_mat1))
  df.melted1 <- melt(data_means1, id = "ncoef_list")
  
  #Plot
  q1 <- ggplot() + #ggplot(data = df.melted, aes(x = ncoef_list, y = value, color = variable)) +
    geom_point(data = df.melted1, aes(x = ncoef_list, y = value, color = variable), inherit.aes = FALSE) +
    geom_line(data = df.melted, aes(x = ncoef_list, y = value, color = variable)) +
    #Tree
    geom_ribbon(data=se.cluster_forest,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#2A9D8F",alpha=0.2, inherit.aes = FALSE) + 
    geom_ribbon(data=se.random_forest,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#F4A261",alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data=se.multi_forest,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#254653",alpha=0.2, inherit.aes = FALSE) +
    #Forest
    geom_ribbon(data=se.cluster_tree,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#af8dc3",alpha=0.2, inherit.aes = FALSE) + 
    geom_ribbon(data=se.random_tree,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#ca0020",alpha=0.2, inherit.aes = FALSE) +
    geom_ribbon(data=se.multi_tree,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#0571b0",alpha=0.2, inherit.aes = FALSE) +
    geom_line(data = df.melted, aes(x = ncoef_list, y = value, color = variable)) +
    theme_classic() + scale_y_continuous(ylab, limits = ylim, breaks  = scales::pretty_breaks(n = 5)) + 
    scale_x_continuous(xlab, limits = xlim, breaks  = scales::pretty_breaks(n = 5)) + labs(title = title) +
    theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)),
          plot.title = element_text(size = rel(1.3)), legend.text=element_text(size=rel(1.3)), 
          legend.title=element_text(size=rel(1.5)), legend.position = "right") +
    scale_color_manual(values=c("#af8dc3","#ca0020", "#0571b0", 
                                "#2A9D8F","#F4A261", "#254653"), 
                       labels = c("Tree: Cluster", "Tree: Random", "Tree: Multi",
                                  "Forest: Cluster", "Forest: Random", "Forest: Multi"), name="Method")
  return(q1)
}
k_list <- c(2, 5, 8, seq(10, 80, 10))

read.csv("~/Desktop/Research 2020/Cluster runs/treeweighting/cs_means_cluster_tree.csv")
treeweighting <- multi.plot.tree(ncoef_list = k_list, filename = "Research 2020/Cluster runs/treeweighting", title = "Weighting Trees vs. Weighting Forests", 
                                 xlab = "Value of k", ylab = "% change in RMSE from Merged", xlim = c(0, 80), ylim = c(-80, 0),n = 100)
show(treeweighting)
