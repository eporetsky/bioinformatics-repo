# R script that enables generating dotplot with your own GO enrichment data
# Find enrichplot information here: https://bioconductor.org/packages/release/bioc/html/enrichplot.html
# The following parser for enrichplot was written by me: Elly Poretsky

library(elliptic)
farey_finder <- function(decimal) {
  # enrichplot requires requires integer-type fraction
  # Find the closest fraction to the decimal using Farey Sequence
  # Maybe doesn't support 1, might need to fix the end loop
  clng <- ceiling(decimal)
  decimal <- decimal/clng
  fr <- farey(50, print=FALSE, give.series = T)
  ln <- dim(fr)[1]
  ln
  for(n in 1:ln){
    temp <- as.vector(fr[n,])
    if(temp[1]/temp[2] > decimal){
      return(paste((clng*temp[1]),temp[2],sep="/"))
    }
    if(n==ln){
      return("1/1")
    }
  }
}

library(enrichplot)
library(DOSE)
# enrichplot is not very customizable so I had to look for ways around
# First, generate the example data provided by the package
data(geneList)
de <- names(geneList)[1:100]
x <- enrichDO(de)

# Next, load your GO enrichment dataset. See provided example
dfg = read.table("enrichplot_dotplot_data", header=T, sep="\t")
# replace the enrichment fold-change number with an integer 
# fraction using the above farey_finder function. 
dfg[,"fc"] <- unlist((lapply(dfg[,"fc"], farey_finder)))
# keep the columns that will be used in in the enrichplot::dotplot function
colnames(x@result)
dfge = dfg[,c("GO", "description", "fc", "n", "p.val","p.adj_bonferroni", "p.adj_bonferroni", "genes", "m")]
# rename the selected columns so it matches the x@result variable
# the data inside some columns doesn't match: no q-values and bg-ratios
# but it is not relevant for the plot I made so I just kept them
colnames(dfge) <- colnames(x@result)
# finally - change x@result to your own GO enrichment data
x@result <- dfge

# Ready to make the enrichplot::dotplot plot
# Further modifications can be done once you can generate the plot
dotplot(x, showCategory=sum(dfg[,"p.adj_bonferroni"]<0.05, na.rm = TRUE))

