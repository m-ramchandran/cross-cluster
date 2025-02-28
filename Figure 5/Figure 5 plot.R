library(rcompanion)
library(ggplot2)
library(reshape2)
library(ggpubr)

#Gene expression
################################
################################

#For pre-clustering strategy:
##########################
#For mutation outcome (continuous)
errors_lgg <- read.csv("~/Desktop/Research 2020/Cluster runs/lgg_gene/errors_total_250.csv")
means_mat <- as.matrix(cbind("Sample Weighted" = errors_lgg$Sample_Weighted, "Stack Ridge" = errors_lgg$Stack_ridge))
apply(errors_lgg, 2, median)

means_mat[,2][which(means_mat[,2] > 400)] <- 0 

df.melted <- melt(means_mat)

q1 <- ggplot(data = df.melted, aes(x = Var2, y = value, group = Var2, fill = Var2)) + geom_boxplot() +#geom_boxplot(fill = fill, colour = line)  + 
  theme_classic() + xlab("Method") + ylab("% change from Merged") + theme(legend.position = "none") + theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3))) +
  scale_fill_brewer(palette = "Set1")    + labs(title = "Performance per iteration") + #ylim(-100, 500) + 
  scale_y_continuous(breaks = c(seq(-100, 0, by = 25), seq(0, 200, by = 100)), limits = c(-100, 100)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black")
show(q1)


## bar plot with means + sds 
#if restricting to the covariates used to cluster over: 
means_lgg <- sapply(1:ncol(means_mat), function(i) median(means_mat[,i][which(-100 < means_mat[,i] & means_mat[,i] < 500)]))
#sds_lgg <- sapply(1:ncol(means_mat), function(i) sd(means_mat[,i][which(-100 < means_mat[,i] & means_mat[,i] < 500)]))


sds_lgg <- t(sapply(1:ncol(means_mat), 
                    function(i) {
                      gM <- groupwiseMedian(data       = as.data.frame(means_mat),
                                            conf       = 0.95,
                                            var = colnames(as.data.frame(means_mat))[i],
                                            R          = 5000,
                                            percentile = TRUE,
                                            bca        = FALSE,
                                            basic      = FALSE,
                                            normal     = FALSE,
                                            wilcox     = FALSE,
                                            digits     = 3)
                      return(c(gM$Percentile.lower, gM$Percentile.upper))
                    }
))


means_lgg <- apply(means_mat,2, median)

plot_means.c <- data.frame(Method = c("Sample Weighted", "Stack Ridge"), 
                           Means = means_lgg,
                           Sds_lower = sds_lgg[,1],
                           Sds_upper = sds_lgg[,2])
plot_means.c$Method <- factor(plot_means.c$Method, levels = plot_means.c$Method)
c1 <- ggplot(plot_means.c, aes(x=Method, y = Means, fill = Method)) + geom_bar(stat = "identity") +
  theme_classic() + geom_errorbar(aes(ymin=Sds_lower, ymax= Sds_upper), width=.2) + scale_fill_brewer(palette = "Set1") + 
  ylab("Median % change")  + theme(legend.position = "none") + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5))) + 
  scale_y_continuous(breaks = seq(-100, 0, by = 20), limits = c(-100, 0)) +
  labs(title = "Median performance") #labs(title = "Performance of approaches on Adult Diffuse Grade Glioma dataset")
show(c1)
annotate_figure(c1, top = text_grob("Ensembling using pre-discovered clusters",size = 20, face = "bold"))

library(ggpubr)
gb <- ggarrange(q1, c1, labels = c("B.1", "B.2"), ncol=1, nrow=2)
show(gb)
gb <- annotate_figure(gb, top = text_grob("Continuous: Mutation count",size = 15, face = "bold"))


###########
#For Grade outcome (Binary)
library("ggsci")

errors_lgg <- rbind(read.csv("~/Desktop/Research 2020/Cluster runs/lgg_gene_binary/errors_total_250.csv"))
means_mat <- as.matrix(cbind("Sample Weighted" = errors_lgg$Sample_Weighted, "Stack Ridge" = errors_lgg$Stack_ridge))
apply(errors_lgg, 2, median)

means_mat[,2][which(means_mat[,2] > 400)] <- 0 

df.melted <- melt(means_mat)

z1 <- ggplot(data = df.melted, aes(x = Var2, y = value, group = Var2, fill = Var2)) + geom_boxplot() +#geom_boxplot(fill = fill, colour = line)  + 
  theme_classic() + xlab("Method") + ylab("% change from Merged") + theme(legend.position = "none") + theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3))) +
  #scale_fill_brewer(palette = "Set2")    + 
  scale_fill_manual(values=c("cadetblue3", "slateblue1")) +
  labs(title = "Performance per iteration") + #ylim(-100, 500) + 
  scale_y_continuous(breaks = c(seq(-40, 0, by = 10), seq(100, 200, by = 100)), limits = c(-40, 0)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black")
show(z1)


## bar plot with means + sds 
#if restricting: 
means_lgg <- sapply(1:ncol(means_mat), function(i) median(means_mat[,i][which(-100 < means_mat[,i] & means_mat[,i] < 500)]))
#sds_lgg <- sapply(1:ncol(means_mat), function(i) sd(means_mat[,i][which(-100 < means_mat[,i] & means_mat[,i] < 500)]))

sds_lgg <- t(sapply(1:ncol(means_mat), 
                    function(i) {
                      gM <- groupwiseMedian(data       = as.data.frame(means_mat),
                                            conf       = 0.95,
                                            var = colnames(as.data.frame(means_mat))[i],
                                            R          = 5000,
                                            percentile = TRUE,
                                            bca        = FALSE,
                                            basic      = FALSE,
                                            normal     = FALSE,
                                            wilcox     = FALSE,
                                            digits     = 3)
                      return(c(gM$Percentile.lower, gM$Percentile.upper))
                    }
))


means_lgg <- apply(means_mat,2, median)

plot_means.c <- data.frame(Method = c("Sample Weighted", "Stack Ridge"), 
                           Means = means_lgg,
                           Sds_lower = sds_lgg[,1],
                           Sds_upper = sds_lgg[,2])
plot_means.c$Method <- factor(plot_means.c$Method, levels = plot_means.c$Method)
d1 <- ggplot(plot_means.c, aes(x=Method, y = Means, fill = Method)) + geom_bar(stat = "identity") +
  theme_classic() + geom_errorbar(aes(ymin=Sds_lower, ymax= Sds_upper), width=.2) + 
  scale_fill_manual(values=c("cadetblue3", "slateblue1")) +
  #scale_fill_brewer(palette = "Set2") + 
  ylab("Median % change")  + theme(legend.position = "none") + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5))) + 
  scale_y_continuous(breaks = seq(-40, 0, by = 10), limits = c(-40, 0)) +
  labs(title = "Median performance") #labs(title = "Performance of approaches on Adult Diffuse Grade Glioma dataset")
show(d1)

library(ggpubr)
zd <- ggarrange(z1, d1, labels = c("C.1", "C.2"), ncol=1, nrow=2)
show(zd)
zd <- annotate_figure(zd, top = text_grob("Binary: Tumor Grade",size = 15, face = "bold"))


final_preset <- ggarrange(gb, zd, ncol = 2, nrow = 1)
final_preset <- annotate_figure(final_preset, top = text_grob("Gene Expression Data" ,size = 20, face = "bold"))
show(final_preset)



##########################
##########################
##########################
#Molecular profiling data

#For mutation outcome - with full merged
errors_lgg <- read.csv("~/Desktop/Research 2020/Cluster runs/lgg_mutation_fullmerged/errors_total_250.csv")

means_mat <- as.matrix(cbind("Sample Weighted" = errors_lgg$Sample_Weighted, "Stack Ridge" = errors_lgg$Stack_ridge, "Subset Merged" = errors_lgg$Merged))

#If removing outliers
#outliers <- boxplot(means_mat, plot=FALSE)$out
#means_mat <- means_mat[-which(means_mat %in% outliers),]

means_mat[,2][which(means_mat[,2] > 400)] <- 0 
df.melted <- melt(means_mat)

e1 <- ggplot(data = df.melted, aes(x = Var2, y = value, group = Var2, fill = Var2)) + geom_boxplot() +#geom_boxplot(fill = fill, colour = line)  + 
  theme_classic() + xlab("Method") + ylab("% change from Merged") + theme(legend.position = "none") + theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3))) +
  scale_fill_brewer(palette = "Set1")    + labs(title = "Performance per iteration") + #ylim(-100, 500) + 
  scale_y_continuous(breaks = c(seq(-100, 100, by = 50), seq(100, 500, by = 100), seq(700, 1200, by = 500)), limits = c(-100, 550)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black")
show(e1)


## bar plot with means + sds 
#if restricting simply to the covariates used to cluster over: 
means_lgg <- sapply(1:ncol(means_mat), function(i) median(means_mat[,i][which(-100 < means_mat[,i] & means_mat[,i] < 500)]))
#sds_lgg <- sapply(1:ncol(means_mat), function(i) sd(means_mat[,i][which(-100 < means_mat[,i] & means_mat[,i] < 500)]))

#if not restricting: 
#means_lgg <- colMeans(means_mat)
#sds_lgg <-  apply(means_mat, 2, sd)

sds_lgg <- t(sapply(1:ncol(means_mat), 
                    function(i) {
                      gM <- groupwiseMedian(data       = as.data.frame(means_mat),
                                            conf       = 0.95,
                                            var = colnames(as.data.frame(means_mat))[i],
                                            R          = 5000,
                                            percentile = TRUE,
                                            bca        = FALSE,
                                            basic      = FALSE,
                                            normal     = FALSE,
                                            wilcox     = FALSE,
                                            digits     = 3)
                      return(c(gM$Percentile.lower, gM$Percentile.upper))
                    }
))


plot_means.c <- data.frame(Method = c("Sample Weighted", "Stack Ridge", "Subset Merged"), 
                           Means = means_lgg,
                           Sds_lower = sds_lgg[,1],
                           Sds_upper = sds_lgg[,2])
plot_means.c$Method <- factor(plot_means.c$Method, levels = plot_means.c$Method)
v1 <- ggplot(plot_means.c, aes(x=Method, y = Means, fill = Method)) + geom_bar(stat = "identity") +
  theme_classic() + geom_errorbar(aes(ymin=Sds_lower, ymax= Sds_upper), width=.2) + scale_fill_brewer(palette = "Set1") + 
  ylab("Median % change")  + theme(legend.position = "none") + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), plot.title = element_text(size = rel(1.5))) + 
  scale_y_continuous(breaks = seq(-40, 0, by = 10), limits = c(-40, 1)) +
  labs(title = "Median performance") #labs(title = "Performance of approaches on Adult Diffuse Grade Glioma dataset")
show(v1)


####################
#Final plot 
ev <- ggarrange(e1, v1, labels = c("A.1", "A.2"), ncol=1, nrow=2)
ev <- annotate_figure(ev, top = text_grob("Continuous: Mutation Count",size = 15, face = "bold"))
ev <- annotate_figure(ev, top = text_grob("Molecular Profiling Data",size = 20, face = "bold"))
show(ev)

final_lgg <- ggarrange(ev, final_preset, ncol = 2, nrow = 1, widths = c(1, 2))
show(final_lgg)