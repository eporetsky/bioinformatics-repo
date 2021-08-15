#!/usr/bin/env Rscript
#args = commandArgs(trailingOnly=TRUE)

# Based on the following repository: https://github.com/eporetsky/parallel_rTASSEL
# Temporarily removed support for parallel GNU implementation as I am working on this file.

# This script is able to calculate the GWAS p-values using rTASSEL and taking the 
# association data p-values, bootstrap it fairly quickly, and re-calculate the 
# SNP marker FDR q-values in an attempt to control low-p-value inflation in an 
# un-biased method.

# Bootstrapping based FDR q-value calculation is currently done using the CIT R package:
# https://www.rdocumentation.org/packages/cit/versions/2.2/topics/cit-package
# Millstein J, Chen GK, Breton CV. 2016. cit: hypothesis testing software for mediation analysis in genomic applications. Bioinformatics. btw135. PMID: 27153715. 
# Millstein J, Zhang B, Zhu J, Schadt EE. 2009. Disentangling molecular relationships with a causal inference test. BMC Genetics, 10:23.

# A single phenotype is replicated n times and randomly sampled. The whole
# phenotype data is run on rTASSEL-FAST mode (which might be better than my
# parallel-rTASSEL implementation) and the p-values of the original observation
# along the n permuations are passed to CIT to calculate q-values.

# According to documentation, CIT provides a conservative FDR q-value adjustment
# if not all p-values are provided, but only a a filtered set of p-values following
# any type of consistent threshold. This helps to avoid JAVA heap memory problems.

# Need sgl to work: sudo apt install libgsl-dev
# devtools::install_github("USCbiostats/cit")
library(cit)
library(rTASSEL)

options(java.parameters = c("-Xmx2g", "-Xms1g"))
rTASSEL::startLogger(fullPath = getwd(), fileName = "RTassel.log")

# Datafrane of the phenotype data
df <- read.table("data.tsv", header=T, sep="\t", row.names=1)
# Number of permutations to run
trait_name <- "trait_name"
col_num <- which(colnames(df)==trait_name)
n_perm <- 1000
# Make a permuation dataframe that is made of n_perm duplicated columns
perm_df <- do.call(cbind, replicate(n_perm, df[,col_num], simplify=FALSE))
# Sample each column to get a random combinations of the values 
# https://stackoverflow.com/questions/6422273/how-to-randomize-or-permute-a-dataframe-rowwise-and-columnwise
perm_df <- as.data.frame(sapply(1:ncol(perm_df), function(col) perm_df[,col]<<-sample(perm_df[,col])))
rownames(perm_df) <- rownames(df)

# Load the 5 PCA Covariates data separately
pc_df <- read.table("5PCs_MDS.txt", header=T, sep="\t", row.names=1, skip=2)
# merge phenotype and covariates by row names
merged <- merge(perm_df, pc_df, by = 0)  
# finally merge the original trait observation data
merged <- cbind(merged, df[,col_num])
# Change the name of the original data to the correct trait name
colnames(merged)[length(merged)] <- colnames(df)[col_num]
# Create the rTASSEL phenotype data from scratch. Notice the attribute vector.
phenoObj <- rTASSEL::readPhenotypeFromDataFrame(phenotypeDF=merged, 
                                                taxaID="Row.names", 
                                                attributeTypes=c(rep(c("data"), n_perm), rep(c("covariate"), 5), "data"))

# Load the rTASSEL genotype/phenotype data with the phenoObj
tasGenoPheno <- rTASSEL::readGenotypePhenotype(
  genoPathOrObj = "genotype.hmp.txt",
  phenoPathDFOrObj = phenoObj
)

# Run the rTASSEL fastAssociation pipeline on all n+1 traits
tasGLM <- rTASSEL::assocModelFitter(
  tasObj = tasGenoPheno,
  formula = . ~ ., # formula that runs all traits at the same time
  fitMarkers = TRUE,
  kinship = NULL,
  fastAssociation = T,
  maxP=0.001
)

# number of observations - number of p-values remained after threshold filtering
nobs <- length(which(tasGLM$FastAssociation$Trait==trait_name))
# ntests is the number of times a p-value calclated, in my case 232977
# ntests is different than nobs because nobs if after p-value filtering
# nobs will also be different for each of the sampled traits.
ntests <- 232977
# Make an empty df for the observed trait p-values. nobs - number of observations
obs = as.data.frame(matrix(NA,nrow=nobs,ncol=2))
names(obs) = c("ID","pvalue")
# ID column is just an ordered list of numbers
obs[,"ID"] = 1:nobs
# Get the list of unique permutation names
perm_names <- unique(tasGLM$FastAssociation$Trait)
# Make a list of vectors that will hold the p value of each permutations
perml = vector('list',n_perm)
count <- 1
for(perm in perm_names){
  if(perm==trait_name){
    nobs <- length(which(tasGLM$FastAssociation$Trait==perm))
    obs["pvalue"] <- tasGLM$FastAssociation[which(tasGLM$FastAssociation$Trait==perm),]$p[1:nobs]
  }
  else{
    nobs <- length(which(tasGLM$FastAssociation$Trait==perm))
    rslts = as.data.frame(matrix(NA,nrow=nobs, ncol=2))
    names(rslts) = c("ID","pvalue")
    rslts[,"ID"] = 1:nobs
    rslts["pvalue"] = tasGLM$FastAssociation[which(tasGLM$FastAssociation$Trait==perm),]$p[1:nobs]
    perml[[count]] <- rslts
    count <- count + 1
  }
}

# Calculate the boostrapped FDR q-values using the CIT function 
q_vals <- fdr.q.perm(obs$pvalue,perml,"pvalue",232977)
# Keep only the SNPs relevant to the original trait observation
tasGLM$FastAssociation <- tasGLM$FastAssociation[which(tasGLM$FastAssociation$Trait == trait_name), ]
# Replace the p-values with the q-values
tasGLM$FastAssociation$p <- q_vals
png(paste("plots/parallelGLM_qvals_",trait_name,".png", sep=""))
manhattanPlot(
  assocStats = tasGLM$FastAssociation,
  trait      = trait_name,
  threshold  = -log10(0.01) # Set a different p-value threshold, not Bonferroni
)
dev.off()