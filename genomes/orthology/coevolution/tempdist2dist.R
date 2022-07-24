library(reshape2)
setwd(paste(getwd(),"/pMT/",sep=""))
fileNames <- Sys.glob("dist/*.tempdist")
for (file in fileNames){
  df <- read.table(file, row.names=1)
  names <- row.names(df)
  names <- gsub("^[^_]*_", "", names)
  colnames(df) <- names
  row.names(df) <- names
  df1 <- as.matrix(df)
  df1 <- df1 + 0.001
  melted <- setNames(melt(df1), c('org1', 'org2', 'distance'))
  write.table(melted,paste(substr(file,1,17),"dist", sep=""), sep="\t", row.names=F, quote=F, col.names=F)
}