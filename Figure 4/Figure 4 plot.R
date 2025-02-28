##############################################################
#Distribution of coefficients plots
##############################################################
library(ggplot2)
library(reshape2)
library(ggpubr)

#Cluster approach
cluster_20 <- read.csv("~/Desktop/Research 2020/Cluster runs/coefs/cs_stack_cluster_20.csv")[,-1]
cluster_80 <- read.csv("~/Desktop/Research 2020/Cluster runs/coefs/cs_stack_cluster_80.csv")[,-1]

list_cluster20 <- c(as.vector(unlist(cluster_20)))
max_20 <- apply(X=cluster_20, MARGIN=1, FUN=max)
clust20 <- c(list_cluster20[list_cluster20 %in% max_20], list_cluster20[!(list_cluster20 %in% max_20)])
max_ind20 <- c(rep("Max weight per ensemble", length(max_20)), rep("Remaining weights", (length(clust20) - length(max_20))))

list_cluster80 <- c(as.vector(unlist(cluster_80)))
max_80 <- apply(X=cluster_80, MARGIN=1, FUN=max)
clust80 <- c(list_cluster80[list_cluster80 %in% max_80], list_cluster80[!(list_cluster80 %in% max_80)])
max_ind80 <- c(rep("Max weight per ensemble", length(max_80)), rep("Remaining weights", (length(clust80) - length(max_80))))

coefs_cluster <- data.frame(k = factor(c(rep("k = 20", each = length(clust20)),
                                         rep("k = 80", each = length(clust80)))), 
                            coefs = c(clust20, clust80),
                            Weights = c(max_ind20, max_ind80))


p1 <- ggplot(coefs_cluster, aes(x=k, y=coefs, fill = Weights)) + geom_boxplot()+  
  xlab("Number of subsections")  + theme_bw() + scale_fill_manual(values = c("#35978f", "#9970ab")) +
  labs(title = "Cluster") + scale_y_continuous(name = "", breaks=seq(0, 4, 1), limits = c(0, 4)) + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), 
        plot.title = element_text(size = rel(1.5)), 
        legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5)))
show(p1)


#Random approach
random_20 <- read.csv("~/Desktop/Research 2020/Cluster runs/coefs/cs_stack_random_20.csv")[,-1]
random_80 <- read.csv("~/Desktop/Research 2020/Cluster runs/coefs/cs_stack_random_80.csv")[,-1]

list_random20 <- c(as.vector(unlist(random_20)))
maxr_20 <- apply(X=random_20, MARGIN=1, FUN=max)
random20 <- c(list_random20[list_random20 %in% maxr_20], list_random20[!(list_random20 %in% maxr_20)])
maxr_ind20 <- c(rep("Max weight per ensemble", length(maxr_20)), rep("Remaining weights", (length(random20) - length(maxr_20))))

list_random80 <- c(as.vector(unlist(random_80)))
maxr_80 <- apply(X=random_80, MARGIN=1, FUN=max)
random80 <- c(list_random80[list_random80 %in% maxr_80], list_random80[!(list_random80 %in% maxr_80)])
maxr_ind80 <- c(rep("Max weight per ensemble", length(maxr_80)), rep("Remaining weights", (length(random80) - length(maxr_80))))

coefs_random <- data.frame(k = factor(c(rep("k = 20", each = length(random20)),
                                        rep("k = 80", each = length(random80)))), 
                           coefs = c(random20, random80),
                           Weights = c(maxr_ind20, maxr_ind80))


p2 <- ggplot(coefs_random, aes(x=k, y=coefs, fill = Weights)) + geom_boxplot()+  
  xlab("Number of subsections") + theme_bw() + scale_fill_manual(values = c("#35978f", "#9970ab"))  +
  labs(title = "Random") + scale_y_continuous(name = "", breaks=seq(0, .7, .1), limits = c(0, .7)) + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), 
        plot.title = element_text(size = rel(1.5)), 
        legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5)))
show(p2)


#Multi approach
multi_20 <- read.csv("~/Desktop/Research 2020/Cluster runs/coefs/cs_stack_multi_20.csv")[,-1]
multi_80 <- read.csv("~/Desktop/Research 2020/Cluster runs/coefs/cs_stack_multi_80.csv")[,-1]

list_multi20 <- c(as.vector(unlist(multi_20)))
maxm_20 <- apply(X=multi_20, MARGIN=1, FUN=max)
multi20 <- c(list_multi20[list_multi20 %in% maxm_20], list_multi20[!(list_multi20 %in% maxm_20)])

maxm_ind20 <- c(rep("Max weight per ensemble", length(maxm_20)), rep("Remaining weights", (length(multi20) - length(maxm_20))))

list_multi80 <- c(as.vector(unlist(multi_80)))
maxm_80 <- apply(X=multi_80, MARGIN=1, FUN=max)
multi80 <- c(list_multi80[list_multi80 %in% maxm_80], list_multi80[!(list_multi80 %in% maxm_80)])
maxm_ind80 <- c(rep("Max weight per ensemble", length(maxm_80)), rep("Remaining weights", (length(multi80) - length(maxm_80))))

coefs_multi <- data.frame(k = (c(rep("5", each = length(multi20)),
                                 rep("5", each = length(multi80)))), 
                          coefs = c(multi20, multi80),
                          Weights = c(maxm_ind20, maxm_ind80))


p3 <- ggplot(coefs_multi, aes(x=k, y=coefs, fill = Weights)) + geom_boxplot()+  
  xlab("Number of true clusters") + theme_bw() + scale_fill_manual(values = c("#35978f", "#9970ab")) +
  labs(title = "Multi") + scale_y_continuous(name = "", breaks=seq(0, .7, .1), limits = c(0, .7))  + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), 
        plot.title = element_text(size = rel(1.5)), 
        legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5)))
show(p3)

library(ggpubr)
#Composite figure
pcomp <- ggarrange(p1, p2, p3, ncol = 3, labels = c("A", "B", "C"), align = "h", nrow=1, common.legend = TRUE, legend="bottom")
annotate_figure(pcomp, left = text_grob("Value of weights", rot = 90, size = 14), 
                top = text_grob("Distribution of Stacking Weights", size = 20, face = "bold"))