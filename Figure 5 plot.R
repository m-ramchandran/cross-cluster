library(ggpubr)
library(ggplot2)
library(reshape2)

##############################################################
#Bias and variance plots
##############################################################

biasvar_100 <- read.csv("~/Desktop/Research 2020/Cluster runs/biasvar/cs_errors_multi100.csv")
biasvar_500 <- read.csv("~/Desktop/Research 2020/Cluster runs/biasvar/cs_errors_multi500.csv")
biasvar_1000 <- read.csv("~/Desktop/Research 2020/Cluster runs/biasvar/cs_errors_multi1000.csv")


biasvar.100 <- as.data.frame(mapply(c, biasvar_100[,1:3], biasvar_100[,4:6]))
biasvar.100 <- cbind(rep(c("Merged", "Multi"), each = 100), biasvar.100)
colnames(biasvar.100) <- c("Method","Bias Squared", "Variance", "MSE")
biasvar.100 <- as.data.frame(biasvar.100)


biasvar.500 <- cbind(rep(c("Merged", "Multi"), each = 100), as.data.frame(mapply(c, biasvar_500[,1:3], biasvar_500[,4:6])))
colnames(biasvar.500) <- c("Method","Bias Squared", "Variance", "MSE")
biasvar.500 <- as.data.frame(biasvar.500)

biasvar.1000 <- cbind(rep(c("Merged", "Multi"), each = 100), as.data.frame(mapply(c, biasvar_1000[,1:3], biasvar_1000[,4:6])))
colnames(biasvar.1000) <- c("Method","Bias Squared", "Variance", "MSE")
biasvar.1000 <- as.data.frame(biasvar.1000)


biasvar <- as.data.frame(rbind(cbind("id" = "n = 500", biasvar.100), cbind("id" = "n = 2500", biasvar.500), cbind("id" = "n = 5000", biasvar.1000)))
biasvar$id <- factor(biasvar$id, levels = levels(as.factor(biasvar$id))[c(2,1,3)], ordered = TRUE)

df.melted <- melt(biasvar, id.vars  = c(1, 2))

b1 <- ggplot(data = df.melted, aes(x=variable, y=value, fill = Method)) + geom_boxplot() + 
  facet_wrap( ~ id, scales="free") +  xlab("Performance metric") + ylab("Value") + theme_bw() + 
  ylim(0, 50) +  labs(title = "Performance per iteration") + scale_fill_brewer(palette = "Paired") + 
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1)), 
        plot.title = element_text(size = rel(1.5)), 
        legend.text=element_text(size=rel(1.3)), legend.title=element_text(size=rel(1.5)))
show(b1)

annotate_figure(b1,top = text_grob("Bias-Variance Decomposition", size = 20, face = "bold"))
