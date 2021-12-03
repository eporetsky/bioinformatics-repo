BiocManager::install("universalmotif")
install.packages("stringi", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
install.packages("stringr", dependencies=TRUE, INSTALL_opts = c('--no-lock')) # might not be required
BiocManager::install("Biostrings", force = TRUE)

library(universalmotif)
library(Biostrings)

# Load all the JASPAR promoter motif data
jaspar <- read_jaspar("JASPAR2022_CORE_plants_non-redundant_pfms_jaspar.txt")
jdf <- summarise_motifs(jaspar)
jdf <- jdf[c("name", "altname")]

# JASPAR has a list of representative motifs for multiple TF families that can be used
family_reps <- read.table("JASPAR_2022_matrix_clustering_plants_CORE_central_motifs_IDs.tab", sep="\t")
family_motifs <- filter_motifs(jaspar, altname=family_reps$V3)

degs <- readDNAStringSet("degs.fasta")
enriched <- enrich_motifs(family_motifs, degs, max.p=0.05, max.q=0.05, max.e=1, nthreads=32) # without ctrl fasta
merged <- merge(enriched, jdf, by.x=c('motif'), by.y=c('name'))
merged <- merged[order(merged$Qval), ]
write.table(merged, "enriched_nCtrl.tsv", sep="\t")

ctrl <- readDNAStringSet("ctrl.fasta")
enriched <-  enrich_motifs(family_motifs, degs, ctrl, max.p=0.05, max.q=0.05, max.e=1, RC=T, nthreads=32) # with ctrl fasta
merged <- merge(enriched, jdf, by.x=c('motif'), by.y=c('name'))
merged <- merged[order(merged$Qval), ]
write.table(merged, "enriched_yCtrl.tsv", sep="\t")

