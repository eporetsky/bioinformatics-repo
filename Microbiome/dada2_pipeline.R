#!/usr/bin/env Rscript

# DADA2 Pipeline for 16S rRNA Amplicon Analysis
# This script processes all samples together for optimal ASV inference

library(dada2)
library(ggplot2)
library(phyloseq)
library(DECIPHER)
library(vegan)      # For PCoA
library(DESeq2)    # For differential abundance

# Read sample metadata
metadata <- read.table("samples.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$accession

# Set paths
path <- "ramdisk/"
output_dir <- "dada2_output"
dir.create(output_dir, showWarnings = FALSE)

# Get list of fastq files
fnFs <- sort(list.files(path, pattern = "_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_2.fastq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

cat("Found", length(sample.names), "samples\n")
cat("Sample names:", paste(sample.names, collapse = ", "), "\n")

# ============================================================================
# 1. Quality Check
# ============================================================================
cat("\n=== Step 1: Quality Profiles ===\n")
# Plot quality profiles for first few samples
pdf(file.path(output_dir, "quality_profiles_forward.pdf"), width = 10, height = 6)
print(plotQualityProfile(fnFs[1:min(4, length(fnFs))]))
dev.off()

pdf(file.path(output_dir, "quality_profiles_reverse.pdf"), width = 10, height = 6)
print(plotQualityProfile(fnRs[1:min(4, length(fnRs))]))
dev.off()

cat("Quality profile plots saved to", output_dir, "\n")

# ============================================================================
# 2. Filter and Trim
# ============================================================================
cat("\n=== Step 2: Filter and Trim ===\n")
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Standard filtering parameters for 16S V3-V4 region (2x250bp or 2x300bp)
# Adjust truncLen based on your quality profiles
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(230, 100),
                     maxN = 0,
                     maxEE = c(2, 2),
                     truncQ = 2,
                     rm.phix = TRUE,
                     compress = TRUE,
                     multithread = TRUE)

cat("Filtering results:\n")
print(head(out))
write.csv(out, file.path(output_dir, "filtering_stats.csv"))

# ============================================================================
# 3. Learn Error Rates
# ============================================================================
cat("\n=== Step 3: Learn Error Rates ===\n")
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Plot error rates
pdf(file.path(output_dir, "error_rates_forward.pdf"), width = 10, height = 10)
print(plotErrors(errF, nominalQ = TRUE))
dev.off()

pdf(file.path(output_dir, "error_rates_reverse.pdf"), width = 10, height = 10)
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

cat("Error rate plots saved to", output_dir, "\n")

# ============================================================================
# 4. Dereplication
# ============================================================================
cat("\n=== Step 4: Dereplication ===\n")
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# ============================================================================
# 5. Sample Inference
# ============================================================================
cat("\n=== Step 5: Sample Inference ===\n")
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

cat("DADA2 inferred sequences for first sample:\n")
print(dadaFs[[1]])

# ============================================================================
# 6. Merge Paired Reads
# ============================================================================
cat("\n=== Step 6: Merge Paired Reads ===\n")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

cat("Merging results for first sample:\n")
print(head(mergers[[1]]))

# ============================================================================
# 7. Construct Sequence Table
# ============================================================================
cat("\n=== Step 7: Construct Sequence Table ===\n")
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
cat("Sequence length distribution:\n")
print(table(nchar(getSequences(seqtab))))

# Save sequence table
saveRDS(seqtab, file.path(output_dir, "seqtab.rds"))

# ============================================================================
# 8. Remove Chimeras
# ============================================================================
cat("\n=== Step 8: Remove Chimeras ===\n")
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", 
                                     multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

cat(sprintf("Chimeras removed: %.1f%% of sequences\n", 
            100 * (1 - sum(seqtab.nochim) / sum(seqtab))))

# Save non-chimeric sequence table
saveRDS(seqtab.nochim, file.path(output_dir, "seqtab_nochim.rds"))

# ============================================================================
# 9. Track Reads Through Pipeline
# ============================================================================
cat("\n=== Step 9: Track Reads ===\n")
getN <- function(x) sum(getUniques(x))
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
cat("Read tracking through pipeline:\n")
print(track)
write.csv(track, file.path(output_dir, "track_reads.csv"))

# ============================================================================
# 10. Assign Taxonomy
# ============================================================================
cat("\n=== Step 10: Assign Taxonomy ===\n")
cat("NOTE: This requires taxonomy database files.\n")
cat("Download Silva or other database using:\n")
cat("wget https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz\n")
cat("wget https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz\n")

# Check if taxonomy files exist
silva_train <- "silva_nr99_v138.1_train_set.fa.gz"
silva_species <- "silva_species_assignment_v138.1.fa.gz"

if (file.exists(silva_train)) {
  cat("Assigning taxonomy using Silva database...\n")
  taxa <- assignTaxonomy(seqtab.nochim, silva_train, multithread = TRUE)
  
  if (file.exists(silva_species)) {
    cat("Adding species-level assignments...\n")
    taxa <- addSpecies(taxa, silva_species)
  }
  
  cat("Taxonomy assignment complete:\n")
  print(head(taxa))
  
  saveRDS(taxa, file.path(output_dir, "taxa.rds"))
  write.csv(taxa, file.path(output_dir, "taxonomy.csv"))
} else {
  cat("WARNING: Taxonomy database not found. Skipping taxonomy assignment.\n")
  cat("Please download Silva database files to the working directory.\n")
  taxa <- NULL
}

# ============================================================================
# 11. Export Results
# ============================================================================
cat("\n=== Step 11: Export Results ===\n")

# Export ASV sequences
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- paste0(">ASV", seq_along(asv_seqs))
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file.path(output_dir, "ASVs.fasta"))

# Export ASV table
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- paste0("ASV", seq_len(nrow(asv_tab)))
write.table(asv_tab, file.path(output_dir, "ASV_table.txt"), 
            sep = "\t", quote = FALSE, col.names = NA)

# Export count table (TSV format for easy import)
write.table(asv_tab, file.path(output_dir, "ASV_counts.tsv"), 
            sep = "\t", quote = FALSE, col.names = NA)

if (!is.null(taxa)) {
  # Export taxonomy with ASV IDs
  taxa_export <- taxa
  row.names(taxa_export) <- paste0("ASV", seq_len(nrow(taxa_export)))
  write.table(taxa_export, file.path(output_dir, "ASV_taxonomy.txt"), 
              sep = "\t", quote = FALSE, col.names = NA)
}

# ============================================================================
# 12. Create Phyloseq Object (if taxonomy available)
# ============================================================================
if (!is.null(taxa)) {
  cat("\n=== Step 12: Create Phyloseq Object ===\n")
  
  # Create sample_data DF (matching order to sequence table)
  valid_samples <- intersect(sample.names, metadata$accession)
  metadata_sub <- metadata[valid_samples, ]

  # Subset sequence table to matched samples only
  sample_idx <- match(valid_samples, sample.names)
  seqtab.nochim <- seqtab.nochim[sample_idx, , drop=FALSE]
  sample.names <- valid_samples

  # Create phyloseq object
  ps <- phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows = FALSE),
    tax_table(taxa),
    sample_data(metadata_sub)
  )
  
  saveRDS(ps, file.path(output_dir, "phyloseq_object.rds"))
  cat("Phyloseq object saved\n")
  
  # Basic summary statistics
  cat("\nPhyloseq object summary:\n")
  print(ps)
  
  # Alpha diversity by group
  if (nsamples(ps) > 1) {
    pdf(file.path(output_dir, "alpha_diversity_by_group.pdf"), width = 8, height = 6)
    p <- plot_richness(ps, measures = c("Shannon", "Simpson", "Observed"), x = "group", color = "group") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
    dev.off()
    cat("Alpha diversity (by group) plot saved\n")
  }
  # PCoA ordination (Bray-Curtis, color by group)
  ord <- ordinate(ps, method = "PCoA", distance = "bray")
  pdf(file.path(output_dir, "PCoA_by_group.pdf"), width = 8, height = 6)
  p_ord <- plot_ordination(ps, ord, color = "group") + geom_point(size = 3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p_ord)
  dev.off()
  cat("PCoA plot saved\n")
  # Differential abundance testing (group vs Col-0 only)
  # Creates a DESeq2 object directly from phyloseq
  if ("Col-0" %in% ps@sam_data$group) {
    dds <- phyloseq_to_deseq2(ps, ~ group)
    dds <- DESeq(dds, fitType = "parametric")
    for (g in setdiff(unique(metadata_sub$group), "Col-0")) {
      res <- results(dds, contrast = c("group", g, "Col-0"))
      write.csv(as.data.frame(res), file = file.path(output_dir, paste0("DESeq2_differential_", g, "_vs_Col-0.csv")))
    }
    cat("Differential abundance testing complete. Results saved to CSV\n")
  }
}

cat("\n=== DADA2 Pipeline Complete! ===\n")
cat("All results saved to:", output_dir, "\n")
cat("\nKey output files:\n")
cat("  - ASVs.fasta: FASTA file of all ASV sequences\n")
cat("  - ASV_counts.tsv: Count table (samples x ASVs)\n")
cat("  - ASV_taxonomy.txt: Taxonomic assignments\n")
cat("  - track_reads.csv: Read counts through each step\n")
cat("  - phyloseq_object.rds: R object for downstream analysis\n")

