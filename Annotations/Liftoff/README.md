# Maize B73 v3 / v4 / v5 gene concordance tables

These tables link gene models across three B73 reference assemblies and annotation releases: **AGPv3 (5b+, `GRMZM*`)** , **Gramene v4 (`Zm00001d*`)** , and **NAM v5 (`Zm00001eb*`)** . Annotations were projected between assemblies with [Liftoff](https://github.com/agshumate/Liftoff) (six hub projections: each version as target once). Genes were linked on each hub by **CDS best-isoform overlap** (≥50% query-relative) on lifted coordinates, with links merged **bidirectionally** (both Liftoff directions) to recover matches where annotation boundaries differ between releases. Liftoff used **full gene models** (including UTR); only the **overlap test** uses CDS intervals. For v3, only **`GRMZM*`** gene models were kept (~36,129 genes); **~3,026 BAC/TE-style models** (`AC*`, `EF*`, etc., ~7.7% of v3 genes) were excluded because they link promiscuously and inflate transitive groups.

### Why CDS overlap (not full gene span)

Many v3 **`GRMZM*`** models have **long or inflated UTRs** that span what v4 and v5 annotate as **separate neighboring genes**. When overlap is measured on the **full gene span** (gene + UTR), one v3 model can overlap several distinct v4/v5 loci after Liftoff, and transitive grouping then **merges unrelated genes into one row**. Measuring overlap on **CDS** (best-matching transcript isoform per query gene) restricts links to **coding loci** and avoids UTR bridges. In a side-by-side run on the same Liftoff output, full-gene linking produced **~12% more multi-copy v5 groups** (4,438 vs 3,919) and **~12% more transitive `multi_mapping` rows** (3,031 vs 2,695) than CDS linking, with little change to strict 1:1 calls. **These published tables use CDS overlap.**

**Files in this directory:** `combined_v5.tsv` — one row per v5 gene (`B73_id` group key); **use as the primary table**. `combined_transitive.tsv` — one row per connected component across all v3/v4/v5 links (union-find); useful for gene families. `hub_v3.tsv`, `hub_v4.tsv`, `hub_v5.tsv` — concordance keyed on each assembly. `multi_group_genes.tsv` — genes appearing on more than one v5-anchored row. Columns: `B73_id`, `5b+`, `Zm00001d`, `Zm00001eb`, `status` (where present). Cell values are gene IDs, comma-separated when ambiguous, or `NONE`.

### Source assemblies and annotations

| Version | Assembly FASTA | Annotation GFF3 |
|---------|----------------|-----------------|
| v3 (5b+) | `B73_RefGen_v3.fa` | `Zea_mays.AGPv3.22.gff3` |
| v4 (Zm00001d) | `Zm-B73-REFERENCE-GRAMENE-4.0.fa` | `Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3` |
| v5 (Zm00001eb) | `Zm-B73-REFERENCE-NAM-5.0.fa` | `Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3` |

Nuclear chromosomes 1–10 only (Mt/Pt and unplaced scaffolds excluded from filtered GFF inputs).
