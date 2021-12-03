BiocManager::install("universalmotif")

library(universalmotif)

# Get the streme (meme suite) enrichment motif analysis results 
streme <- read_meme("streme.txt", skip = 0, readsites = FALSE, readsites.meta = FALSE)
# Download the Jaspar database from their website and read it
jaspar <- read_jaspar("JASPAR2022_CORE_plants_non-redundant_pfms_jaspar.txt")

# Extract the TF family names from the Jaspar file
jdf <- summarise_motifs(jaspar)
jdf <- jdf[c("name", "altname")]

# Compare all the jaspar motif to the enriched streme motifs
comps <- compare_motifs(c(streme, jaspar), 1:length(streme), method="PCC", normalise.scores=TRUE,
                        min.mean.ic=.0, min.overlap=5, score.strat="a.mean")

# Merge the comparison results with the Jaspar TF names
merged <- merge(comps, jdf, by.x=c('target'), by.y=c('name'))
merged <- merged[c("subject","altname","target","score","logPval","Eval")]
merged[order(merged$logPval), ]
write.table(merged, "streme_results.tsv", sep="\t")