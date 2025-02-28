####################
#Multistudy
ncoef_list <-c(2,5, 8, 10, 15, seq(20, 80, 10))
means_cluster <- means_multi <- means_random <- means_RF <- means_ridge <- means_nnet <- as.data.frame(array(0, c(length(ncoef_list), 13)))

colnames_total <- c("Merged", "Unweighted", "CS_Weighted",
                    "Stack_noint", "Stack_noint_norm", "Stack_int",
                    "SS_noint", "SS_noint_norm", "SS_int", "Stack_lasso", "SS_lasso", "Stack_ridge", "SS_ridge")

filename1 <- "250reps/multistudy"


means_cluster <- read.csv(paste0("~/Desktop/", filename1, "/cs_means_cluster.csv"))[-c(13, 14), ]
means_multi <- read.csv(paste0("~/Desktop/", filename1, "/cs_means_multi.csv"))[-c(13, 14), ]
means_random <- read.csv(paste0("~/Desktop/", filename1, "/cs_means_random.csv"))[-c(13, 14), ]
means_nnet <- read.csv(paste0("~/Desktop/", filename1, "/cs_means_nnet.csv"))[-c(13, 14), ]
means_nnet_cluster <- c(3.969782, 3.906399, 3.861816, 3.832086, 3.814205, 3.723049, 3.581760, 
                        3.550759, 3.491533, 3.428731, 3.375482, 3.326065)

means_RF_merged <- do.call("rbind", replicate(nrow(means_cluster), mean(means_cluster$Merged), simplify = FALSE))
means_RF_multi <- do.call("rbind", replicate(nrow(means_multi), mean(means_multi$Stack_ridge), simplify = FALSE))
means_nnet_merged <- do.call("rbind", replicate(nrow(means_nnet), mean(means_nnet$Merged), simplify = FALSE))
means_nnet_multi <- do.call("rbind", replicate(nrow(means_nnet), mean(means_nnet$Stack_ridge), simplify = FALSE))

means_mat <- as.matrix(cbind(means_RF_merged, means_cluster$Stack_ridge, means_RF_multi, means_random$Stack_ridge,
                             means_nnet_merged, means_nnet_cluster, means_nnet_multi))

data_means <- data.frame(cbind(ncoef_list, means_mat))
df.melted <- melt(data_means, id = "ncoef_list")

means_mat1 <- as.matrix(cbind(means_RF_merged, means_RF_multi,
                              means_nnet_merged, means_nnet_multi))
data_means1 <- data.frame(cbind(ncoef_list, means_mat1))
df.melted1 <- melt(data_means1, id = "ncoef_list")



q1 <- ggplot(data = df.melted, aes(x = ncoef_list, y = value, color = variable)) +
  geom_point()  + geom_line() +
  theme_classic() + ylab("Average RMSE") + 
  scale_x_continuous("Value of k", limits = c(0, 80), breaks  = scales::pretty_breaks(n = 8)) + 
  labs(title = "Multistudy setting: Comparing performance of various learners") +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)), 
        plot.title = element_text(size = rel(1.5)), 
        legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5))) +
  scale_color_manual(values=c("blue","deepskyblue2", "seagreen1", "black",
                              "#ef8a62","#762a83", "#b2182b"), 
                     labels = c("RF:Merged", "RF: Cluster", "RF: Multi", "RF:Random", 
                                "Neural Net: Merged", "Neural Net: Cluster", "Neural Net: Multi"), name="Method")
show(q1)