library(ggplot2)

# Data with columns: Treatment, Genotype and Value
dat <- read.table("data.csv", sep=",", header=T)

# Assuming there is only 2 unique factors in the data
dat$Value.f <- factor(ifelse(dat$Value == "label1", 1, 0), 
                labels = c("label1", "label2"))

plt <- ggplot(dat, aes(x=Treatment, y="Log2(FC)", fill=Genotype)) +
  geom_boxplot() + 
  scale_x_discrete(name = "") +
    ggtitle("") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          text = element_text(size = 12, family = "Tahoma", face = "bold"),
          axis.title = element_text(size=12, face="bold"),
          axis.text.x=element_text(size = 12, face="bold")) +
    facet_grid(.~Value.f) + 
  scale_color_manual(values=c("blue", "red"))