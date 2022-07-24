#BiocManager::install("edgeR")

library(edgeR)

# Load the raw count data
counts <- read.table("raw.counts", header=TRUE, row.names=1, sep="\t")
# Prepare the relevant group variable. Same order as count file
groups <- c(rep("Ctrl", 4), rep("Treatment1", 4), rep("Treatment2", 4))
ds <- DGEList(counts=counts, group=factor(groups))
ds$samples$lib.size <- colSums(ds$counts)
ds <- calcNormFactors(ds)
ds <- estimateCommonDisp(ds)
comparisons <- c("2H", "4H", "8H", "24H", "48H")
for(comp in comparisons){
  DE <- exactTest(ds, pair=c("Ctrl",comp))$table
  write.table(DE, paste("FC/sorghum_timecourse",comp,".tsv", sep=""), sep="\t", quote=F)
}
#  ----------------------------------------------------

#  ----------------------------------------------------Sorghum global data
count <- read.table("files/Counts/Tetreault_sorghum_counts.txt", header=TRUE, row.names=1, sep="\t")
# Maize res_sus data
groups <- c(rep("RTx2783_D5_Ctrl", 3), rep("RTx2783_D5_Aphid", 3), rep("RTx2783_D10_Ctrl", 3), rep("RTx2783_D10_Aphid", 3), rep("RTx2783_D15_Ctrl", 3), rep("RTx2783_D15_Aphid", 3),
            rep("BCK60_D5_Ctrl", 3), rep("BCK60_D5_Aphid", 3), rep("BCK60_D10_Ctrl", 3), rep("BCK60_D10_Aphid", 3), rep("BCK60_D15_Ctrl", 3), rep("BCK60_D15_Aphid", 3))
ds <- DGEList(counts=count, group=factor(groups))
ds$samples$lib.size <- colSums(ds$counts)
ds <- calcNormFactors(ds)
ds <- estimateCommonDisp(ds)
comparisons <- c("2H", "4H", "8H", "24H", "48H")
for(comp in comparisons){
  DE <- exactTest(ds, pair=c("Ctrl",comp))$table
  write.table(DE, paste("FC/sorghum_timecourse",comp,".tsv", sep=""), sep="\t", quote=F)
}
#  ----------------------------------------------------

#  ----------------------------------------------------Our sorghum field data
count <- read.table("files/Counts/sorghum_field_counts.txt", header=TRUE, row.names=1, sep="\t")
groups <- c(rep("F", 4), rep("WT", 4))
ds <- DGEList(counts=count, group=factor(groups))
ds$samples$lib.size <- colSums(ds$counts)
ds <- calcNormFactors(ds)
ds <- estimateCommonDisp(ds)
comparisons <- c("2H", "4H", "8H", "24H", "48H")
for(comp in comparisons){
  DE <- exactTest(ds, pair=c("Ctrl",comp))$table
  write.table(DE, paste("FC/sorghum_timecourse",comp,".tsv", sep=""), sep="\t", quote=F)
}
#  ----------------------------------------------------

#  ----------------------------------------------------Maize aphid data
count <- read.table("Tzin_maize_timecourse_counts.txt", header=TRUE, row.names=1, sep="\t")
groups <- c(rep("Ctrl1", 5), rep("2H", 5), rep("4H", 5), rep("8H", 5),
            rep("24H", 5), rep("48H", 5), rep("96H", 5), rep("Ctrl2", 5))
ds <- DGEList(counts=count, group=factor(groups))
ds$samples$lib.size <- colSums(ds$counts)
ds <- calcNormFactors(ds)
ds <- estimateCommonDisp(ds)
comparisons <- c("2H", "4H", "8H", "24H", "48H", "96H")
for(comp in comparisons){
  DE <- exactTest(ds, pair=c("Ctrl2",comp))$table
  write.table(DE, paste("FC/Tzin_maize_timecourse_",comp,".tsv", sep=""), sep="\t", quote=F)
}
#  ----------------------------------------------------

#  ----------------------------------------------------Maize res_sus data
count <- read.table("Pingault_maize_counts.txt", header=TRUE, row.names=1, sep="\t")
groups <- c(rep("Mp708-0h-L", 3), rep("Mp708-24h-L", 3), rep("Mp708-0h-R", 3), rep("Mp708-24h-R", 3),
            rep("Tx601-0h-L", 3), rep("Tx601-24h-L", 3), rep("Tx601-0h-R", 3), rep("Tx601-24h-R", 3))
ds <- DGEList(counts=count, group=factor(groups))
ds$samples$lib.size <- colSums(ds$counts)
ds <- calcNormFactors(ds)
ds <- estimateCommonDisp(ds)
DE <- exactTest(ds, pair=c("Mp708-0h-L","Mp708-24h-L"))$table
write.table(DE, paste("FC/Mp708_Res_FC.tsv", sep=""), sep="\t", quote=F)
DE <- exactTest(ds, pair=c("Tx601-0h-L","Tx601-24h-L"))$table
write.table(DE, paste("FC/Tx601_Sus_FC.tsv", sep=""), sep="\t", quote=F)
#  ----------------------------------------------------

#  ----------------------------------------------------Sorghum res_sus data
count <- read.table("Kiani_sorghum_counts.txt", header=TRUE, row.names=1, sep="\t")
groups <- c(rep("37-07-ctrl_2wk", 3), rep("37-07-aphid_2wk", 3), rep("37-07-ctrl_6wk", 3), rep("37-07-aphid_6wk", 3),
            rep("44-20-ctrl_2wk", 3), rep("44-20-aphid_2wk", 3), rep("44-20-ctrl_6wk", 3), rep("44-20-aphid_6wk", 3))
ds <- DGEList(counts=count, group=factor(groups))
ds$samples$lib.size <- colSums(ds$counts)
ds <- calcNormFactors(ds)
ds <- estimateCommonDisp(ds)
DE <- exactTest(ds, pair=c("37-07-ctrl_2wk","37-07-aphid_2wk"))$table
write.table(DE, paste("FC/37-07_Res_FC.tsv", sep=""), sep="\t", quote=F)
DE <- exactTest(ds, pair=c("44-20-ctrl_2wk","44-20-aphid_2wk"))$table
write.table(DE, paste("FC/44-20_Sus_FC.tsv", sep=""), sep="\t", quote=F)

DE <- exactTest(ds, pair=c("37-07-ctrl_6wk","37-07-aphid_6wk"))$table
write.table(DE, paste("FC/37-07_Res_6wks_FC.tsv", sep=""), sep="\t", quote=F)
DE <- exactTest(ds, pair=c("44-20-ctrl_6wk","44-20-aphid_6wk"))$table
write.table(DE, paste("FC/44-20_Sus_6_wks_FC.tsv", sep=""), sep="\t", quote=F)

DE <- exactTest(ds, pair=c("44-20-aphid_2wk","37-07-aphid_2wk"))$table
write.table(DE, paste("FC/44-20_vs_37-07_2wks_FC.tsv", sep=""), sep="\t", quote=F)
DE <- exactTest(ds, pair=c("44-20-aphid_6wk","37-07-aphid_6wk"))$table
write.table(DE, paste("FC/44-20_vs_37-07_6wks_FC.tsv", sep=""), sep="\t", quote=F)

#  ----------------------------------------------------

#  ----------------------------------------------------maize aphid time-course data
count <- read.table("Tzin_maize_timecourse_counts.txt", header=TRUE, row.names=1, sep="\t")
groups <- c(rep("Ctrl1", 5), rep("2H", 5), rep("4H", 5), rep("8H", 5),
            rep("24H", 5), rep("48H", 5), rep("96H", 5), rep("Ctrl2", 5))
ds <- DGEList(counts=count, group=factor(groups))
ds$samples$lib.size <- colSums(ds$counts)
ds <- calcNormFactors(ds)
ds <- estimateCommonDisp(ds)
comparisons <- c("2H", "4H", "8H", "24H", "48H")
for(comp in comparisons){
  DE <- exactTest(ds, pair=c("Ctrl",comp))$table
  write.table(DE, paste("FC/sorghum_timecourse",comp,".tsv", sep=""), sep="\t", quote=F)
}
#  ----------------------------------------------------

#  ----------------------------------------------------

# Calculate the exact test  ----------------------------------------------------
DE2 <- exactTest(ds, pair=c("Ctrl","2H"))$table
DE2["FDR"] <- p.adjust(DE2$PValue, method="BH")
DE4 <- exactTest(ds, pair=c("Ctrl","4H"))$table
DE4["FDR"] <- p.adjust(DE4$PValue, method="BH")
DE8 <- exactTest(ds, pair=c("Ctrl","8H"))$table
DE8["FDR"] <- p.adjust(DE8$PValue, method="BH")
DE24 <- exactTest(ds, pair=c("Ctrl","24H"))$table
DE24["FDR"] <- p.adjust(DE24$PValue, method="BH")
DE48 <- exactTest(ds, pair=c("Ctrl","48H"))$table
DE48["FDR"] <- p.adjust(DE48$PValue, method="BH")
#DE96 <- exactTest(ds, pair=c("Ctrl2","96H"))$table
#DE96["FDR"] <- p.adjust(DE96$PValue, method="BH")

#summary(decideTestsDGE(exactTest(dSb, pair=c("Ctrl","2H")), adjust.method="BH", p.value=0.05, lfc=1))

fdr_thresh <- Inf
fc_trhesh <- 0
DE2 <- subset(DE2[DE2["FDR"]<fdr_thresh & abs(DE2["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
DE4 <- subset(DE4[DE4["FDR"]<fdr_thresh & abs(DE4["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
DE8 <- subset(DE8[DE8["FDR"]<fdr_thresh & abs(DE8["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
DE24 <- subset(DE24[DE24["FDR"]<fdr_thresh & abs(DE24["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
DE48 <- subset(DE48[DE48["FDR"]<fdr_thresh & abs(DE48["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
#DE96 <- subset(DE96[DE96["FDR"]<fdr_thresh & abs(DE96["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))

# Keep log2FC

fdr_thresh <- Inf
fc_trhesh <- 0
DE2 <- subset(DE2[DE2["FDR"]<fdr_thresh & abs(DE2["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
DE4 <- subset(DE4[DE4["FDR"]<fdr_thresh & abs(DE4["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
DE8 <- subset(DE8[DE8["FDR"]<fdr_thresh & abs(DE8["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
DE24 <- subset(DE24[DE24["FDR"]<fdr_thresh & abs(DE24["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
DE48 <- subset(DE48[DE48["FDR"]<fdr_thresh & abs(DE48["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
#DE96 <- subset(DE96[DE96["FDR"]<fdr_thresh & abs(DE96["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))


# Reduce takes a binary function and a list of data items and successively 
# applies the function to the list elements in a recursive fashion.
DEcomb <- Reduce(
  function(x, y) {
    temp <- merge(x, y, all = TRUE, by="row.names") %>%
            column_to_rownames(., var = "Row.names")
    
  }, list(DE2, DE4, DE8, DE24, DE48)#, DE96)
)
colnames(DEcomb) <- c("2H", "4H", "8H", "24H", "48H")#, "96H")
DEcomb[is.na(DEcomb)] <- 0
write.table(DEcomb, "FC/sorghum_timecourse__DEcombAll_pvals.tsv", sep="\t", quote=F)
write.table(DEcomb, "FC/sorghum_timecourse__DEcombAll_log2FC.tsv", sep="\t", quote=F)

#########################################################################################

countZm <- read.table("coco_zm.csv", header=TRUE, row.names=1, sep=",")
ZmDataGroups <- c(rep("Ctrl", 5), rep("2H", 5), rep("4H", 5), rep("8H", 5), rep("24H", 5), rep("48H", 5), rep("96H", 5))
dZm <- DGEList(counts=countZm, group=factor(ZmDataGroups))
dZm$samples$lib.size <- colSums(dZm$counts)
dZm <- calcNormFactors(dZm)
plotMDS(dZm, method="logFC", col=as.numeric(dZm$samples$group))
legend("bottomleft", as.character(unique(dZm$samples$group)), col=1:6, pch=20)
dZm <- estimateCommonDisp(dZm)

# Calculate the exact test  ----------------------------------------------------
write.table(exactTest(dZm, pair=c("Ctrl","2H")), "Final_FC_Results/Zm2H-edgeRexact.tsv", sep="\t", quote=F) # compare groups 1 and 2
write.table(exactTest(dZm, pair=c("Ctrl","4H")), "Final_FC_Results/Zm4H-edgeRexact.tsv", sep="\t", quote=F) # compare groups 1 and 2
write.table(exactTest(dZm, pair=c("Ctrl","8H")), "Final_FC_Results/Zm8H-edgeRexact.tsv", sep="\t", quote=F) # compare groups 1 and 2
write.table(exactTest(dZm, pair=c("Ctrl","24H")), "Final_FC_Results/Zm24H-edgeRexact.tsv", sep="\t", quote=F) # compare groups 1 and 2
write.table(exactTest(dZm, pair=c("Ctrl","48H")), "Final_FC_Results/Zm48H-edgeRexact.tsv", sep="\t", quote=F) # compare groups 1 and 2
write.table(exactTest(dZm, pair=c("Ctrl","96H")), "Final_FC_Results/Zm96H-edgeRexact.tsv", sep="\t", quote=F) # compare groups 1 and 2
write.table(exactTest(dZm, pair=c("2H","48H")), "Final_FC_Results/Zm2H-48H-edgeRexact.tsv", sep="\t", quote=F) # compare groups 1 and 2


DE2 <- exactTest(dZm, pair=c("Ctrl","2H"))$table
DE2["FDR"] <- p.adjust(DE2$PValue, method="BH")
write.table(DE2, "Final_FC_Results/Zm2H-edgeRexact.tsv", sep="\t", quote=F)
DE4 <- exactTest(dZm, pair=c("Ctrl","4H"))$table
DE4["FDR"] <- p.adjust(DE4$PValue, method="BH")
write.table(DE4, "Final_FC_Results/Zm4H-edgeRexact.tsv", sep="\t", quote=F)
DE8 <- exactTest(dZm, pair=c("Ctrl","8H"))$table
DE8["FDR"] <- p.adjust(DE8$PValue, method="BH")
write.table(DE8, "Final_FC_Results/Zm8H-edgeRexact.tsv", sep="\t", quote=F)
DE24 <- exactTest(dZm, pair=c("Ctrl","24H"))$table
DE24["FDR"] <- p.adjust(DE24$PValue, method="BH")
write.table(DE24, "Final_FC_Results/Zm24H-edgeRexact.tsv", sep="\t", quote=F)
DE48 <- exactTest(dZm, pair=c("Ctrl","48H"))$table
DE48["FDR"] <- p.adjust(DE48$PValue, method="BH")
write.table(DE48, "Final_FC_Results/Zm48H-edgeRexact.tsv", sep="\t", quote=F)
DE96 <- exactTest(dZm, pair=c("Ctrl","96H"))$table
DE96["FDR"] <- p.adjust(DE96$PValue, method="BH")
write.table(DE96, "Final_FC_Results/Zm96H-edgeRexact.tsv", sep="\t", quote=F)

summary(decideTestsDGE(exactTest(dZm, pair=c("Ctrl","48H")), adjust.method="BH", p.value=0.05, lfc=1))

fdr_thresh <- Inf
fc_trhesh <- 0
DE2 <- subset(DE2[DE2["FDR"]<fdr_thresh & abs(DE2["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
DE4 <- subset(DE4[DE4["FDR"]<fdr_thresh & abs(DE4["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
DE8 <- subset(DE8[DE8["FDR"]<fdr_thresh & abs(DE8["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
DE24 <- subset(DE24[DE24["FDR"]<fdr_thresh & abs(DE24["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
DE48 <- subset(DE48[DE48["FDR"]<fdr_thresh & abs(DE48["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))
DE96 <- subset(DE96[DE96["FDR"]<fdr_thresh & abs(DE96["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, FDR))

# Reduce takes a binary function and a list of data items and successively 
# applies the function to the list elements in a recursive fashion.
ZmDEcomb <- Reduce(
  function(x, y) {
    temp <- merge(x, y, all = TRUE, by="row.names") %>%
      column_to_rownames(., var = "Row.names")
    
  }, list(DE2, DE4, DE8, DE24, DE48, DE96)
)
colnames(ZmDEcomb) <- c("2H", "4H", "8H", "24H", "48H", "96H")
ZmDEcomb[is.na(ZmDEcomb)] <- 0
write.table(ZmDEcomb, "Final_FC_Results/ZmDEcombAll_filtered_log2FC.tsv", sep="\t", quote=F)

# Need to re-run all the DE2-48

DE2 <- subset(DE2[DE2["FDR"]<fdr_thresh & abs(DE2["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
DE4 <- subset(DE4[DE4["FDR"]<fdr_thresh & abs(DE4["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
DE8 <- subset(DE8[DE8["FDR"]<fdr_thresh & abs(DE8["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
DE24 <- subset(DE24[DE24["FDR"]<fdr_thresh & abs(DE24["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
DE48 <- subset(DE48[DE48["FDR"]<fdr_thresh & abs(DE48["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))
DE96 <- subset(DE96[DE96["FDR"]<fdr_thresh & abs(DE96["logFC"])>=fc_trhesh,], select=-c(logCPM, PValue, logFC))


# Reduce takes a binary function and a list of data items and successively 
# applies the function to the list elements in a recursive fashion.
ZmDEcomb <- Reduce(
  function(x, y) {
    temp <- merge(x, y, all = TRUE, by="row.names") %>%
      column_to_rownames(., var = "Row.names")
    
  }, list(DE2, DE4, DE8, DE24, DE48, DE96)
)
colnames(ZmDEcomb) <- c("2H", "4H", "8H", "24H", "48H", "96H")
ZmDEcomb[is.na(ZmDEcomb)] <- 0
write.table(ZmDEcomb, "Final_FC_Results/ZmDEcombAll_filtered_FDR.tsv", sep="\t", quote=F)

DE2 <- exactTest(dZm, pair=c("Ctrl","2H"))$table
DE4 <- exactTest(dZm, pair=c("Ctrl","4H"))$table
DE8 <- exactTest(dZm, pair=c("Ctrl","8H"))$table
DE24 <- exactTest(dZm, pair=c("Ctrl","24H"))$table
DE48 <- exactTest(dZm, pair=c("Ctrl","48H"))$table

DE2 <- subset(DE2[DE2["PValue"]<0.05 & abs(DE2["logFC"])>1,], select=-c(logCPM, PValue))
DE4 <- subset(DE4[DE4["PValue"]<0.05 & abs(DE4["logFC"])>1,], select=-c(logCPM, PValue))
DE8 <- subset(DE8[DE8["PValue"]<0.05 & abs(DE8["logFC"])>1,], select=-c(logCPM, PValue))
DE24 <- subset(DE24[DE24["PValue"]<0.05 & abs(DE24["logFC"])>1,], select=-c(logCPM, PValue))
DE48 <- subset(DE48[DE48["PValue"]<0.05 & abs(DE48["logFC"])>1,], select=-c(logCPM, PValue))

# Reduce takes a binary function and a list of data items and successively 
# applies the function to the list elements in a recursive fashion.
ZmDEcomb <- Reduce(
  function(x, y) {
    temp <- merge(x, y, all = TRUE, by="row.names") %>%
      column_to_rownames(., var = "Row.names")
    
  }, list(DE2, DE4, DE8, DE24, DE48)
)
colnames(ZmDEcomb) <- c("2H", "4H", "8H", "24H", "48H")
ZmDEcomb[is.na(ZmDEcomb)] <- 0
# First I need to write the Zm files and rename the rows to the Sorghum IDs
write.table(ZmDEcomb, "ZmDEcomb.tsv", sep="\t", quote=F)
ZmDEcomb <- read.table("ZmDEcombConvert2Sobic.tsv.txt", head=T, row.names=1, sep="\t")

ZmDEcomb
SbZmMerged <- merge(SbDEcomb, ZmDEcomb, all = FALSE, by="row.names")

# Calculate the GLM test  ----------------------------------------------------
# GLM estimates of dispersion
design.mat <- model.matrix(~ 0 + dSb$samples$group)
colnames(design.mat) <- levels(dSb$samples$group)
dSb2 <- estimateGLMCommonDisp(dSb,design.mat)
dSb2 <- estimateGLMTrendedDisp(dSb2, design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
dSb2 <- estimateGLMTagwiseDisp(dSb2,design.mat)
#plotBCV(dSb2)
dSb$samples$group


fit <- glmFit(dSb2, design.mat)
DE2 <- glmLRT(fit, contrast=c(0,1,0,0,0,-1))$table
DE4 <- glmLRT(fit, contrast=c(0,0,0,1,0,-1))$table
DE8 <- glmLRT(fit, contrast=c(0,0,0,0,1,-1))$table
DE24 <- glmLRT(fit, contrast=c(1,0,0,0,0,-1))$table
DE48 <- glmLRT(fit, contrast=c(0,0,1,0,0,-1))$table
DElate <- glmLRT(fit, contrast=c(1,0,1,0,0,-1))$table
DEearly <- glmLRT(fit, contrast=c(0,1,0,1,0,-1))$table

DE2 <- subset(DE2[DE2["PValue"]<0.05 & abs(DE2["logFC"])>1,], select=-c(logCPM, LR, PValue))
DE4 <- subset(DE4[DE4["PValue"]<0.05 & abs(DE4["logFC"])>1,], select=-c(logCPM, LR, PValue))
DE8 <- subset(DE8[DE8["PValue"]<0.05 & abs(DE8["logFC"])>1,], select=-c(logCPM, LR, PValue))
DE24 <- subset(DE24[DE24["PValue"]<0.05 & abs(DE24["logFC"])>1,], select=-c(logCPM, LR, PValue))
DE48 <- subset(DE48[DE48["PValue"]<0.05 & abs(DE48["logFC"])>1,], select=-c(logCPM, LR, PValue))
DElate <- subset(DElate[DElate["PValue"]<0.05 & abs(DElate["logFC"])>1,], select=-c(logCPM, LR, PValue))
DEearly <- subset(DEearly[DEearly["PValue"]<0.05 & abs(DEearly["logFC"])>1,], select=-c(logCPM, LR, PValue))

SbDEglmComb <- Reduce(
  function(x, y) {
    temp <- merge(x, y, all = TRUE, by="row.names") %>%
      column_to_rownames(., var = "Row.names")
    
  }, list(DE2, DE4, DE8, DE24, DE48, DElate, DEearly)
)
colnames(SbDEglmComb) <- c("2H", "4H", "8H", "24H", "48H", "Late", "Early")
SbDEglmComb[is.na(SbDEglmComb)] <- 0
# First I need to write the Zm files and rename the rows to the Sorghum IDs
write.table(SbDEglmComb, "SbDEglmComb.tsv", sep="\t", quote=F)
ZmDEcomb <- read.table("ZmDEcombConvert2Sobic.tsv.txt", head=T, row.names=1, sep="\t")







ZmDataGroups
colnames(countZm)
countZm <- read.table("coco_zm.csv", header=TRUE, row.names=1, sep=",")
ZmDataGroups <- c(rep("Ctrl", 5), rep("2H", 5), rep("4H", 5), rep("8H", 5), rep("24H", 5), rep("48H", 5), rep("96H", 5))
dZm <- DGEList(counts=countZm, group=factor(ZmDataGroups))
dZm$samples$lib.size <- colSums(dZm$counts)
dZm <- calcNormFactors(dZm)
plotMDS(dZm, method="logFC", col=as.numeric(dZm$samples$group))
legend("bottomleft", as.character(unique(dZm$samples$group)), col=1:6, pch=20)
dZm <- estimateCommonDisp(dZm)
# GLM estimates of dispersion
design.mat <- model.matrix(~ 0 + dZm$samples$group)
colnames(design.mat) <- levels(dZm$samples$group)
dZm2 <- estimateGLMCommonDisp(dZm,design.mat)
dZm2 <- estimateGLMTrendedDisp(dZm2, design.mat, method="power")
# You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
# The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
dSb <- estimateGLMTagwiseDisp(dZm2,design.mat)
plotBCV(dZm2)



# Calculate average CPM and correlation  ----------------------------------------------------

?read.table
write.table(cpm(dSb), "coco_SbCPM.csv", sep=",", quote=F)
write.table(cpm(dZm), "coco_ZmCPM.csv", sep=",", quote=F)

# Calculate the correlation between each gene pair ------------------------------------------

SbOrth <- read.table("coco_SbCPM_mean_ZmMatchingOrder.csv", header=TRUE, row.names=1, sep=",")
ZmOrth <- read.table("coco_ZmCPM_mean_SbMatchingOrder.csv", header=TRUE, row.names=1, sep=",")
# pvals <- sapply(seq.int(dim(SbOrth)[1]), function(i) cor.test(t(log2(SbOrth[i,])), t(log2(ZmOrth[i,])))$p.value)

test<- cor.test(t(SbOrth[2,]), t(ZmOrth[2,]))
test$estimate
A <- as.matrix(df1)
B <- as.matrix(df2)
pvals <- sapply(seq.int(dim(SbOrth)[1]), function(i) cor.test(t(SbOrth[i,]), t(ZmOrth[i,]), method="spearman")$p.value)
pvals <- p.adjust(pvals, method="BH")

cvals <- sapply(seq.int(dim(SbOrth)[1]), function(i) cor.test(t(SbOrth[i,]), t(ZmOrth[i,]), method="spearman")$estimate)
#cvals <- 1 / cvals
write.table(dfCors, "CoCoResults.csv", sep=",", quote=F)
dfCors <- data.frame("cor" = cvals, "p"=pvals)
?cor.test
p1 <- ggplot(dfCors, aes(-log10(p), cor)) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab("Inverse PCC") + 
  ylab("FDR p-val")
p1


head(rowCors)
p.adjust




countdata <- subset(countdata, select = -c(X24H_3))


# Assign condition (first four are controls, second four contain the expansion)
treatment <- factor(c(rep("Control", 4), 
                      rep("Aphid", 20)
))

time <- factor(c(rep("48H", 4), 
                 rep("2H", 4),
                 rep("4H", 4),
                 rep("8H", 4),
                 rep("24H", 4),
                 rep("48H", 4)
))

coldata <- data.frame(row.names=colnames(countdata), treatment, time)
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, 
                              design = ~ treatment + time)# + treatment:time)
dds <- DESeq(dds, test="LRT", reduced = ~ treatment) 
resultsNames(dds)
write.table(results(dds, contrast=c("time", "48H", "Ctrl")), "LRT3_48H-FC.csv", quote = F)


?resultsNames

dds <- varianceStabilizingTransformation(dds)
coldata
# https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)



dds <- varianceStabilizingTransformation(dds)

BiocManager::install("fission")
library("fission")
data("fission")
ddsTC <- DESeqDataSet(fission, ~ strain + minute + strain:minute)
ddsTC@design

fission@colData

#coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
reducedDesign <- as.formula(~ Ctrl )
reducedDesign <- as.formula(~ 1 )
dds@design
~condition

ddsTC <- DESeqDataSet(fission, ~ treatment + time + treatment:minute)


dds <- DESeq(dds, test="LRT", reduced = reducedDesign)
resultsNames(dds)
?resultsNames
write.table(results(dds, contrast=c("condition", "48H", "Ctrl")), "LRT3_48H-FC.csv", quote = F)
?DESeq
#normalized_counts <- counts(dds, normalized=TRUE)
#write.table(normalized_counts, file="200203_SbTimeCourse_normalizedCounts.csv", sep=",", quote=F, col.names=NA)
#rld <- rlog(dds)

sampleDists <- dist(t(vst(assay(dds, blind=FALSE))))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(countdata)
colnames(sampleDistMatrix) <- colnames(countdata)
library("gplots")
library("RColorBrewer")
colours = colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
hclust.method <- function(x) hclust(x, method="ward.D2")
heatmap.2(sampleDistMatrix, trace="none", col=colours, hclustfun=hclust.method)

pr.hc.s <- hclust(sampleDists, method = "single")
pr.hc.c <- hclust(sampleDists, method = "complete")
pr.hc.a <- hclust(sampleDists, method = "average")
pr.hc.w <- hclust(sampleDists, method = "ward.D2")
# plot them
op <- par(mar = c(0, 4, 4, 2), mfrow = c(2, 2))
dev.off()
plot(pr.hc.s, labels = colnames(countdata), main = "Single", xlab = "")
plot(pr.hc.c, labels = colnames(countdata), main = "Complete", xlab = "")
plot(pr.hc.a, labels = colnames(countdata), main = "Average", xlab = "")
plot(pr.hc.w, labels = colnames(countdata), main = "Ward", xlab = "")

# Analysis of FC with DESeq2 ----------------------------------------------------
# rld <- rlog(dds)
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd)

# Analysis of FC with DESeq2 ----------------------------------------------------
vsdf <- assay(vsd)
drop_cols <- colnames(vsdf)
avg_df <- DataFrame(Ctrl = as.vector(rowMeans(vsdf[,drop_cols[1:4]])),
          H2   = as.vector(rowMeans(vsdf[,drop_cols[5:8]])),
          H4   = as.vector(rowMeans(vsdf[,drop_cols[9:12]])),
          H8   = as.vector(rowMeans(vsdf[,drop_cols[13:16]])),
          H24  = as.vector(rowMeans(vsdf[,drop_cols[17:19]])),
          H48  = as.vector(rowMeans(vsdf[,drop_cols[20:23]]))
)

rownames(avg_df) <- rownames(countdata)


# draw a figure showing time varying gene expression
# in each cluster, overlain with the each clusters 
# mean time series
alpha.long %>%
  ggplot(aes(time, expression, group=gene)) +
  geom_line(alpha=0.25) + 
  geom_line(aes(time, mean.exp, group=NULL,color=cluster),
            data = cluster.means,
            size=1.1) +
  ylim(-2.5, 2.5) +
  facet_wrap(~cluster, ncol=4)



# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds)
res <- results(dds)

write.table(results(dds, contrast=c("condition", "Ctrl", "4H")), "4H-FC.csv", quote = F)
res

# Analysis time-course with DESeq2 ----------------------------------------------
countdata <- read.table("200203_SbTimeCourseNoCtrl.csv", header=TRUE, row.names=1, sep=",")

# Assign condition (first four are controls, second four contain the expansion)
condition <- factor(c(rep("2H", 4),
                      rep("4H", 4),
                      rep("8H", 4),
                      rep("24H", 4),
                      rep("48H", 4)
))
coldata <- data.frame(row.names=colnames(countdata), condition)

# https://hbctraining.github.io/DGE_workshop_salmon/lessons/02_DGE_count_normalization.html
dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

dds

vsd <- vst(dds)
head(assay(vsd), 3)
df <- as.data.frame(colData(dds)[,c("condition")])
vsd

library(RColorBrewer)
(mycols <- brewer.pal(6, "Dark2")[1:length(unique(condition))])
library(gplots)
# Sample distance heatmap
#sampleDists <- as.matrix(dist(t(assay(rld))))
rld <- rlogTransformation(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:1000]
#sampleDists <- assay(vsd)[select,]
sampleDists <- rld[select,]
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], #RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()




BiocManager::install("EBSeq")
library(EBSeq)
EB_norm <- MedianNorm(countdata, alternative = FALSE)
?MedianNorm
head(countdata)
# Run the DESeq pipeline
#dds <- DESeq(dds)
?DESeq
dds<-DESeq(ddfs, test="Wald", fitType="parametric", sfType="poscounts", minReplicatesForReplace=3)
?DESeq
# Plot dispersions
# png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
# dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
rld <- dds
head(assay(rld))
hist(assay(rld))

?varianceStabilizingTransformation
rld <- varianceStabilizingTransformation(dds)
# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(6, "Dark2")[1:length(unique(condition))])
library(gplots)
# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
#png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
#dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
## I like mine better:
rld_pca <- function(rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
#png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
#dev.off()


# Get differential expression results
res <- results(dds)
head(res)
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Write results
write.csv(resdata, file="diffexpr-results.csv")

## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
#png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
#dev.off()