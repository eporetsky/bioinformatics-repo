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
  # Top N ranks by median abundance
  topN <- abund %>% group_by(!!sym(taxrank)) %>% summarize(Total = median(RelAbund)) %>% top_n(N, Total) %>% pull(!!sym(taxrank))
  abund <- abund %>% mutate(Taxon = ifelse((!!sym(taxrank)) %in% topN, (!!sym(taxrank)), "Other"))
  # Combine all 'Other' for each sample (sum abundances)
  abund <- abund %>%
    mutate(Taxon = ifelse((!!sym(taxrank)) %in% topN, (!!sym(taxrank)), "Other")) %>%
    group_by(Sample, Taxon, group) %>%
    summarize(
      Abundance = sum(Abundance),
      RelAbund = sum(RelAbund),
      .groups = "drop"
    )

  # Set proper Taxon levels and palette
  abund$Taxon <- factor(abund$Taxon, levels = c("Other", setdiff(unique(abund$Taxon), "Other")))
  n_col <- length(unique(abund$Taxon))
  palette <- colorRampPalette(brewer.pal(min(9, n_col), "Set1"))(n_col)

  # Order samples by group and sample, using meta
  meta_with_names <- meta %>% tibble::rownames_to_column("Sample")
  if ("group" %in% colnames(meta_with_names)) {
    sample_order <- meta_with_names %>% arrange(group, Sample) %>% pull(Sample)
  } else {
    sample_order <- meta_with_names %>% arrange(Sample) %>% pull(Sample)
  }
  abund$Sample <- factor(abund$Sample, levels = sample_order)

  # By-group faceted plot only
  p <- ggplot(abund, aes(x = Sample, y = RelAbund, fill = Taxon)) +
    geom_bar(stat = "identity", color = "black", size = 0.2) +
    scale_fill_manual(values = palette, drop = FALSE) +
    facet_grid(. ~ group, scales = "free_x", space = "free") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.2, "lines"),
      strip.text.x = element_text(size=10),
      legend.position = "right",
      legend.title = element_text(size=9),
      legend.text = element_text(size=8)
    ) +
    guides(fill = guide_legend(ncol = 1, title = taxrank, reverse = FALSE))

  ggsave(filename, p, width = 12, height = 6)
}

# 1. Stacked Barplots (Family, Top 10 and Top 20)
generate_stacked_barplot(ps, "Family", N = 15, pastel=TRUE, filename="dada2_output/relabund_family_samples_top15.png")
