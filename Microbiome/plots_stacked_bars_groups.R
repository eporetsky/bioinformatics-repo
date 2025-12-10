#!/usr/bin/env Rscript

library(phyloseq)
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(tibble)

# Load Phyloseq object
ps <- readRDS("dada2_output/phyloseq_object.rds")
meta <- data.frame(sample_data(ps))

# A function for relative abundance barplots at a given tax level
generate_stacked_barplot <- function(ps, taxrank, N = 10, pastel = TRUE, filename) {
  agg_ps <- tax_glom(ps, taxrank)
  abund <- psmelt(agg_ps)
  abund <- abund %>% 
    group_by(Sample, !!sym(taxrank)) %>% 
    summarize(Abundance = sum(Abundance)) %>%
    group_by(Sample) %>% mutate(RelAbund = Abundance / sum(Abundance)) %>%
    ungroup()
  abund <- abund %>% left_join(meta %>% tibble::rownames_to_column("Sample"), by = "Sample")
  # Use only Family column throughout
  # Compute top N families
  topN <- abund %>% group_by(Family) %>% summarize(Total = median(RelAbund)) %>%
    top_n(N, Total) %>% pull(Family)

  # Assign topN or 'Other' as Taxon
  abund <- abund %>%
    mutate(Taxon = ifelse(Family %in% topN, Family, "Other"))

  # Sum for each group & taxon
  grouped_abund <- abund %>%
    group_by(group, Taxon) %>%
    summarize(
      Abundance = sum(Abundance),
      .groups = "drop"
    ) %>%
    group_by(group) %>%
    mutate(RelAbund = Abundance / sum(Abundance)) %>%
    ungroup()

  # Set Taxon factor and palette
  grouped_abund$Taxon <- factor(grouped_abund$Taxon, levels = c("Other", setdiff(unique(grouped_abund$Taxon), "Other")))
  n_col <- length(unique(grouped_abund$Taxon))
  palette <- colorRampPalette(brewer.pal(min(9, n_col), "Set1"))(n_col)

  # By-group barplot (one bar per group)
  p <- ggplot(grouped_abund, aes(x = group, y = RelAbund, fill = Taxon)) +
    geom_bar(stat = "identity", color = "black", size = 0.2) +
    scale_fill_manual(values = palette, drop = FALSE) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      legend.title = element_text(size=9),
      legend.text = element_text(size=8)
    ) +
    guides(fill = guide_legend(ncol = 1, title = "Family", reverse = FALSE))

  ggsave(filename, p, width = 8, height = 6)
}

# 1. Stacked Barplots (Family, Top 10 and Top 20)
generate_stacked_barplot(ps, "Family", N = 15, pastel=TRUE, filename="dada2_output/relabund_family_groups_top15.png")
