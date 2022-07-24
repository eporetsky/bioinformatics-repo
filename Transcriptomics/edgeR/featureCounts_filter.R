library(dplyr) # I don't think that this is used

# Read the count file produced by featureCounts 
df <- read.table("counts.txt", sep="\t", header=T, row.names=1, skip=1)

# Remove the featureCounts columns containing the information about genes leaving raw counts
df <- select(df, -c(colnames(df)[1:5]))

# Write table to file
write.table(df, "counts_raw.txt", quote=F, sep="\t")
