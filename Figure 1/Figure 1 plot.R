#####################
#Plots 
library(ggplot2)
library(ggpubr)

#INSTRUCTIONS:
#First run all code in Figure1a_code.R, Figure1b_code.R, and Figure1c_code.R
#This code assumes particular directories for each of the results, but can be easily changed. 

# Plotting function
ratio_plot <- function(filename, plot_title, 
                       ncoef_list = c(2, 4, 8, 16, 32)){
  errors_cluster <- do.call(cbind, lapply(ncoef_list, function(ncoef){
    read.csv(paste0("~/Desktop/Research/Research 2021/rf_theory/", filename,"/cs_errors_cluster", ncoef+1, ".csv"))
  }))
  errors_merged <- do.call(cbind, lapply(ncoef_list, function(ncoef){
    read.csv(paste0("~/Desktop/Research/Research 2021/rf_theory/", filename,"/cs_errors_merged", ncoef+1, ".csv"))
  }))
  
  empirical_ratio <- errors_cluster/errors_merged
   n = 150
  
  plot_data <- cbind.data.frame("k" = ncoef_list, 
                                "mean_ratio" = colMeans(empirical_ratio), 
                                "sd_ratio" = apply(empirical_ratio, 2, sd))
  
  
  # Create the plot
  plot <- ggplot(plot_data, aes(x = k, y = mean_ratio)) +
    geom_point(size = 3) +  # plot the means as points
    geom_ribbon(aes(ymin = (mean_ratio - 1.96*sd_ratio/sqrt(n)), ymax = (mean_ratio + 1.96*sd_ratio/sqrt(n))), 
                fill = "lightblue", alpha = 0.5) +  # add error bars
    geom_line() + 
    geom_hline(yintercept = 1/sqrt(2), linetype = "dashed", color = "blue") + 
    labs(title = plot_title,
         x = "Number of clusters (k)", 
         y = "") +
    theme_minimal() + ylim(.65, .92)
  
  # Display the plot
  return(plot)
}

ratio_plot("uniform", "Uniform clusters")
ratio_plot("gaussian_restrict_update", "Gaussian clusters")
ratio_plot("laplace", "Laplace clusters")

#composite plot
rf_theory_composite <- ggarrange(ratio_plot("uniform", "Uniform clusters"), 
                                 ratio_plot("gaussian_restrict_update", "Gaussian clusters"), 
                                 ratio_plot("laplace", "Laplace clusters"), 
                                 labels = c("A", "B", "C"),
                                 ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
show(rf_theory_composite)

annotate_figure(rf_theory_composite, left = text_grob("Average RMSE ratio", rot = 90, size = 14),
                top = text_grob("Average Ratio of Ensemble RMSE to Merged RMSE", size = 20, face = "bold"))

