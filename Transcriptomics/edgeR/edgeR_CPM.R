#BiocManager::install("edgeR")

library(edgeR)

file_name <- "name"

# Load the raw count data
counts <- read.table(paste(file_name, "_counts.tsv", sep=""), header=TRUE, row.names=1, sep="\t")
# Load the data and process with edgeR
ds <- DGEList(counts=counts)
ds$samples$lib.size <- colSums(ds$counts)
ds <- calcNormFactors(ds)
ds <- estimateCommonDisp(ds)

# Write CPM data to TSV file
write.table(cpm(ds), paste(file_name, "_cpm.tsv", sep=""), sep="\t", quote=F)

# Make an MDS plot of the count data
png(paste(file_name, "_MDS.png", sep=""))
plotMDS(ds, method="logFC")
dev.off()
