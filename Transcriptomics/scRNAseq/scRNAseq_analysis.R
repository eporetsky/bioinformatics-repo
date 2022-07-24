# An attempt to run through the steps of scRNA-seq analysis
# Using conda instead of regular methods because of multiple failures.
# conda install -c conda-forge r-seurat
# conda install -c conda-forge r-stringr
# conda install -c conda-forge r-umap


library(Matrix)
# Convert raw expression data to a sparse matrix
expression <- read.delim("scRNAseq_expression_matrix.csv", header=T, row.names=1, sep=",")
mat_sparse <- Matrix(as.matrix(genes), sparse=TRUE)
writeMM(mat_sparse, "matrix.mtx")

library(dplyr)
library(patchwork)

# Read10X didn't work with my partial gene name file so I saved all names and duplicated column
write.table(data.frame(row.names(genes)), "genes.tsv", sep="\t", quote=F)

library(Seurat)
### I am following this tutorial, step-by-step analysis: 
### https://danhdtruong.com/Seurat-Guided-Clustering-Tutorial/

# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "Seurat/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "scRoots", min.cells = 3, min.features = 200)

# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# Visualize QC metrics as ridge plots
RidgePlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol =1)


# I am skipping the QC and normalization steps because I assume the GEO data is post-processed
# library(scater) # BiocManager::install("scater")
# qc.nCount_RNA <- isOutlier(pbmc$nCount_RNA, log=TRUE, type="both")
# qc.nFeature_RNA  <- isOutlier(pbmc$nFeature_RNA, log=TRUE, type="both")
# pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# pbmc <- NormalizeData(pbmc)

library(ggplot2)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc) + theme(legend.text = element_text(size = 6))
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + theme(legend.text = element_text(size = 6))
plot1 + plot2

pbmc <- ScaleData(pbmc, features = VariableFeatures(object = pbmc)) #VariableFeatures is used to call the highly variable genes from the object.

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc), verbose = FALSE)
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(pbmc, reduction = "pca")
ElbowPlot(pbmc, ndims = 50)

# Here we use the first 50 PCs to construct the neighbor graph.
pbmc <- FindNeighbors(pbmc, dims = 1:50)

# Next, we can apply algorithms to identify “communities” of cells. 
# There are two main community detection algorithm, Louvain and the improved version Leiden. 
# We use the default Louvain while controlling the resolution parameter to adjust the number of clusters.
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc <- RunTSNE(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
plot1 <- DimPlot(pbmc, reduction = "umap")
plot2 <- DimPlot(pbmc, reduction = "tsne")
plot3 <- DimPlot(pbmc, reduction = "pca")
plot1 +  plot2 + plot3 + plot_layout(nrow = 2)
