library(ggplot2)
library(reshape2)
library(ggpubr)

multi.plot <- function(ncoef_list, filename, title, xlab, ylab,xlim, ylim, n){
  paste0("~/Desktop/", filename, "/cs_means_cluster.csv")
  
  colnames_total <- c("Merged", "Unweighted", "CS_Weighted",
                      "Stack_noint", "Stack_noint_norm", "Stack_int",
                      "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")
  
  means_cluster <- read.csv(paste0("~/Desktop/", filename, "/cs_means_cluster.csv"), skip = 1, header = F)[-c(18, 19), -1]
  
  sds_cluster <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_cluster.csv"), skip = 1, header = F)[-c(18, 19), -1]
  
  means_multi <- read.csv(paste0("~/Desktop/", filename, "/cs_means_multi.csv"), header = F, skip = 1)[-c(18, 19), -1]
  sds_multi <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_multi.csv"), skip = 1, header = F)[-c(18, 19), -1]
  
  means_multi <- do.call("rbind", replicate(nrow(means_multi), means_multi[5,], simplify = FALSE))
  sds_multi <- do.call("rbind", replicate(nrow(sds_multi), sds_multi[5,], simplify = FALSE))
  
  
  means_random <- read.csv(paste0("~/Desktop/", filename, "/cs_means_random.csv"), skip = 1, header = F)[-c(18, 19), -1]
  sds_random <- read.csv(paste0("~/Desktop/", filename, "/cs_sds_random.csv"), skip = 1, header = F)[-c(18, 19), -1]
  
  colnames(means_cluster) <- colnames(means_multi) <- colnames(means_random) <- colnames_total
  colnames(sds_cluster) <- colnames(sds_multi) <- colnames(sds_random) <- colnames_total
  
  
  #Creating the plotting data frames
  means_mat <- as.matrix(cbind(means_cluster$Stack_ridge, means_random$Stack_ridge, means_multi$Stack_ridge))
  data_means <- data.frame(cbind(ncoef_list, means_mat))
  df.melted <- melt(data_means, id = "ncoef_list")
  se.cluster <- data.frame(cbind(ncoef_list, means_cluster$Stack_ridge - 1.96*sds_cluster$Stack_ridge/(sqrt(n)),
                                 means_cluster$Stack_ridge + 1.96*sds_cluster$Stack_ridge/(sqrt(n))))
  se.multi <- data.frame(cbind(ncoef_list, means_multi$Stack_ridge - 1.96*sds_multi$Stack_ridge/(sqrt(n)),
                               means_multi$Stack_ridge + 1.96*sds_multi$Stack_ridge/(sqrt(n))))
  se.random <- data.frame(cbind(ncoef_list, means_random$Stack_ridge - 1.96*sds_random$Stack_ridge/(sqrt(n)),
                                means_random$Stack_ridge + 1.96*sds_random$Stack_ridge/(sqrt(n))))
  
  means_mat1 <- as.matrix(cbind(means_cluster$Stack_ridge, means_random$Stack_ridge))
  data_means1 <- data.frame(cbind(ncoef_list, means_mat1))
  df.melted1 <- melt(data_means1, id = "ncoef_list")
  
  #Plot
  q1 <- ggplot() + #ggplot(data = df.melted, aes(x = ncoef_list, y = value, color = variable)) +
    geom_point(data = df.melted1, aes(x = ncoef_list, y = value, color = variable), inherit.aes = FALSE) +
    geom_line(data = df.melted, aes(x = ncoef_list, y = value, color = variable)) +
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

#Binary plots
k_list <- c(2:10, seq(10, 80, 10))
binary_gaussian <- multi.plot(ncoef_list = k_list, filename = "250reps/binarized/linear", title = "Gaussian clusters", 
                              xlab = "Value of k", ylab = "% change in RMSE from Merged", xlim = c(0, 80), ylim = c(-21, 10),n = 250)
show(binary_gaussian)

binary_monte <- multi.plot(ncoef_list = k_list, filename = "250reps/binarized/monte", title = "Non-gaussian clusters", 
                           xlab = "Value of k", ylab = "% change in RMSE from Merged",xlim = c(0, 80), ylim = c(-21, 10), n = 250)
show(binary_monte)


binary_monte <- annotate_figure(binary_monte,
                                top = text_grob("Binary", size = 17, face = "bold"))

#Linear plots

linear_gaussian <- multi.plot(ncoef_list = k_list, filename = "250reps/linear/gaussian", title = "Gaussian clusters", 
                              xlab = "Value of k", ylab = "% change in RMSE from Merged", xlim = c(0, 80), ylim = c(-48, -4),n = 250)
show(linear_gaussian)

linear_monte <- multi.plot(ncoef_list = k_list, filename = "250reps/linear/monte", title = "Non-gaussian clusters", 
                           xlab = "Value of k", ylab = "% change in RMSE from Merged",xlim = c(0, 80), ylim = c(-48, -4), n = 250)
show(linear_monte)


linear_monte <- annotate_figure(linear_monte, top = text_grob("Linear", size = 17, face = "bold"))

#Quadratic plots

quadratic_gaussian <- multi.plot(ncoef_list = k_list, filename = "250reps/quadratic/linear", title = "Gaussian clusters", 
                                 xlab = "Value of k", ylab = "% change in RMSE from Merged", xlim = c(0, 80), ylim = c(-35, 6),n = 500)
show(quadratic_gaussian)

quadratic_monte <- multi.plot(ncoef_list = k_list, filename = "250reps/quadratic/monte", title = "Non-gaussian clusters", 
                              xlab = "Value of k", ylab = "% change in RMSE from Merged",xlim = c(0, 80), ylim = c(-35, 6), n = 250)
show(quadratic_monte)

quadratic_monte <- annotate_figure(quadratic_monte, top = text_grob("Quadratic", size = 17, face = "bold"))



setk <- ggarrange(linear_monte, binary_monte, quadratic_monte, 
                  linear_gaussian, binary_gaussian, quadratic_gaussian,
                  labels = c("A.1", "B.1", "C.1", "A.2", "B.2", "C.2"),
                  ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
show(setk)

annotate_figure(setk,left = text_grob("% change in average RMSE from Merged", rot = 90, size = 14),
                top = text_grob("Varying outcome models", size = 20, face = "bold"))
