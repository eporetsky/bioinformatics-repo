# **WIP:** kmersGWAS Pipeline Setup and Usage

This repository contains scripts and environment setup to help with running the k-mers-based GWAS analysis using the [kmersGWAS package](https://github.com/voichek/kmersGWAS/tree/master/src).

## 1. Environment Setup

1. **Create the conda environment:**
   ```bash
   conda env create -f environment.yml
   conda activate kmersgwas
   ```
2. **Install the R package `mvnpermute` (not available via conda):**
   ```bash
   R
   > install.packages("mvnpermute")
   # Then quit R
   ```

## 2. Download and Set Up kmersGWAS

Run the setup script to download kmersGWAS, KMC, and GEMMA:
```bash
bash 01_download_kmersGWAS.sh
```
This will download all the repository into the current directory with all required binaries.

## 3. Prepare Your Data
- Place all your filtered FASTQ files in a directory (e.g., `fastq/`).
- Ensure your phenotype table is available (e.g., `trait.txt`).

## 4. Run the Pipeline

### Step 2: KMC Counting for Each Individual
This script will loop over all FASTQ files in the specified directory, group them by sample, and process each sample:
```bash
bash 02_kmc_per_individual_template.sh fastq 31 2
```
- Arguments: `<KMER_LENGTH> <THREADS>`
- Output: For each sample, a directory with `kmers_with_strand` and logs.

### Step 3: Combine and Filter k-mers
Edit `kmers_list_paths.txt` to list all `kmers_with_strand` files and sample names (tab-separated). Then run:
```bash
bash 03_combine_and_filter_kmers.sh 31 5 0.2
```
- Arguments: `<KMER_LENGTH> <MAC> <PERCENT>`

### Step 4: Build the k-mers Table
```bash
bash 04_build_kmers_table.sh 31
```
- Arguments: `<KMER_LENGTH>`

### Step 5: Advanced Analysis (Kinship, PLINK, GWAS)
```bash
bash 05_advanced_analysis_templates.sh 31 trait.txt 8
```
- Arguments: `<KMER_LENGTH> <PHENO_FILE> <THREADS>`

## Notes
- Adjust script parameters as needed for your dataset.
- Phenotype file is tab-separated with “accession_id” and “phenotype_value” headers
- Cite the original work: Voichek, Y., & Weigel, D. (2020). Identifying genetic variants underlying phenotypic variation in plants without complete genomes. Nature genetics, 52(5), 534-540.