library(ape)
library(phylogram)

# input csv distance matrix
file_name <- "matrix.dist"

# Read the distance matrix
df_all <- read.table(paste(file_name, ".csv", sep=""), sep=",", header=T, row.names=1)

# Generate a neighbor-joining tree from the distance matrix
nj_all <- nj(as.matrix(df_all))

# Should root the tree in mid-point, not sure if needed
rt_all <- remidpoint(nj_all)

# Convert the tree object to a phylogram dendrogram object
dn_all <- as.dendrogram(rt_all)

# Use phylogram to write the dendrogram object
write.dendrogram(dn_all, file=paste(file_name, ".nwk", sep=""), edges = TRUE)

# Plots the dendrogram
plot(dn_all)