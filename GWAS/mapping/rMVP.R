################## Using rMVP for GWAS mapping ################## 
# For full instructions see: https://github.com/xiaolei-lab/rMVP/
# A fully imputed hapmap file is necessary to run rMVP
# To generate imputed hapmap I used LinkImpute. Find jupyter notebook:
# https://github.com/eporetsky/bioinformatics-assortment/blob/master/GWAS/editors/hapmapLinkImpute.ipynb

devtools::install_github("xiaolei-lab/rMVP")

library(rMVP)

# Full-featured function (Recommended)
MVP.Data(fileHMP="hapmap.txt",
         filePhe="Phenotype.txt",
         sep.phe="\t",
         fileKin=FALSE,
         filePC="PCs.txt",
         #priority="memory",
         #maxLine=10000,
         out="mvp.hmp"
)

# Only convert genotypes
# the genotype data should be fully imputed before using this function
MVP.Data.Hapmap2MVP("hapmap.txt", out='mvp') 

genotype <- attach.big.matrix("mvp.hmp.geno.desc") 
phenotype <- read.table("mvp.hmp.phe",head=TRUE)
map <- read.table("mvp.hmp.geno.map" , head = TRUE)

imMVP <- MVP(
  phe=phenotype,
  geno=genotype,
  map=map,
  nPC.GLM=5,   ##if you have added PCs into covariates, please keep there closed.
  nPC.MLM=3,  ##if you don't want to add PCs as covariates, please comment out the parameters instead of setting the nPC to 0.
  nPC.FarmCPU=3,
  priority="speed",   ##for Kinship construction
  vc.method="BRENT",  ##only works for MLM, BRENT, HE, EMMA
  maxLoop=10,
  method.bin="static",   ## "FaST-LMM", "static" (#only works for FarmCPU)
  threshold=0.05,
  method=c("GLM", "MLM", "FarmCPU")
)