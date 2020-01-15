load(url("http://dmirman.github.io/FAex.Rdata"))

#https://david-d.ncifcrf.gov/helps/functional_annotation.html
#https://www.slideshare.net/framancuso/david-5451863
#http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td.html

library(reshape2)
library(ggplot2)
library(scales)

go_df <- read.csv("GO_terms.csv", header=T)
go_df$description <- factor(go_df$description, levels = go_df$description[order(go_df$DEG_item)])
go_df$type_new <- factor(go_df$type,levels=c("Upregulated","Downregulated"))

g<-ggplot(go_df, aes(description, DEG_item, fill=Fold_enrichment)) + 
  facet_grid(rows=vars(type_new), 
             scales = "free_y", space="free_y",drop=T) + #place the factors in separate facets
  scale_y_continuous(expand=c(0,0)) +
  geom_col() + #make the bars
  coord_flip() + #flip the axes so the test names can be horizontal  
  #define the fill color gradient: blue=positive, red=negative
  scale_fill_gradient2(name = "Fold\nEnrichment", 
        high = "red", 
        limits=c(0,12),
        guide=guide_colorbar(frame.colour = "black", ticks.colour = "black", frame.linewidth = 1,
                             ticks = F, label.position	="left", bin=50)) +
  ylab("Number of genes") +
  xlab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent",colour = NA),
        
        legend.title.align=0.5,
        legend.justification = 'center',
        legend.text=element_text(colour = "black", family = "Helvetica"),
        
        panel.border = element_rect(color = "black", fill = NA, size = 1), 
        strip.background = element_rect(color = "black", size = 1),
        axis.line = element_line(colour = "black"),
        strip.text.y = element_text(angle = -90),
        strip.placement = "outside",
        strip.text = element_text(size=9, lineheight=0.5,colour="black",family="Helvetica",face="bold"),
        axis.text.y = element_text(colour = "black", family = "Helvetica"),
        axis.text.x = element_text(colour = "black", family = "Helvetica")
        )

ggsave("GO_terms.png",g, height = 5, width = 6, dpi=320, bg="white")
