#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Based on the following repository: https://github.com/eporetsky/parallel_rTASSEL
# This script does the same as above and also calculates the lambdaGC stratification
# value. My goal is too find an un-biased measurement for filtering traits that have
# an inflation of low p-values that could get in the way of downstream analysis. 

# Used this as a reference for calculating lambdaGC: 
# https://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
# This has additional ideas for testing if p-values follow a uniform distribution: 
# https://stats.stackexchange.com/questions/110755/calculate-inflation-observed-and-expected-p-values-from-uniform-distribution-in

library(rTASSEL)

# input parameter of number of traits in the phenotype data
col_num <- as.integer(args[1])

options(java.parameters = c("-Xmx2g", "-Xms1g"))
rTASSEL::startLogger(fullPath = getwd(), fileName = "RTassel.log")

tasGenoPheno <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = "genotype.hmp.txt",
  phenoPathDFOrObj = "phenotype.txt"
)

feature_names <- names(rTASSEL::getPhenotypeDF(tasObj = tasGenoPheno))
# Skip first column and last 5 PC1-5 columns to get a vector of trait names
feature_names <- feature_names[2:(length(feature_names)-5)]

formula <- as.formula(paste(feature_names[col_num], "~ ."))

tasGLM <- rTASSEL::assocModelFitter(
  tasObj = tasGenoPheno,
  formula = formula,
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = F,
  maxP=1 # to calculate lambda all SNPs are needed
)

chisq <- qchisq(1-tasGLM$GLM_Stats$p,1)
chisq <- chisq[!is.na(chisq)]
lambdagc <- round(median(chisq)/qchisq(0.5,1), 2)

write.table(tasGLM$GLM_Stats[order(tasGLM$GLM_Stats$p), ], paste("results/parallelGLM_", lambdagc, "_", feature_names[col_num], ".tsv", sep=""), sep="\t", quote=FALSE)

png(paste("plots/parallelGLM_", lambdagc, "_", feature_names[col_num], ".png", sep=""))
manhattanPlot(
  assocStats = tasGLM$GLM_Stats,
  trait      = feature_names[col_num],
  threshold  = -log10(2.2e-7)
)
dev.off()

tasGLM$GLM_Stats <- tasGLM$GLM_Stats[which(tasGLM$GLM_Stats$p < 2.2e-7),] # calculated bonferroni threshold for 233037 SNPs

write.table(tasGLM$GLM_Stats[order(tasGLM$GLM_Stats$p), ], paste("results/parallelGLM_", lambdagc, "_", feature_names[col_num], ".tsv", sep=""), sep="\t", quote=FALSE)
