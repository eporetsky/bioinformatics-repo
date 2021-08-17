library(caret)

df <- read.table("data.csv", header=T, row.names=1, sep=",")
df <- log2(df+1)
df <- df[-caret::nearZeroVar(t(df)),]

write.table(df, "data_nearZeroVar.csv", quote=FALSE, sep=",")
