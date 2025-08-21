#!/bin/bash
set -e

# Usage: ./05_advanced_analysis_templates.sh <KMER_LENGTH> <PHENO_FILE> <THREADS>
# Example: ./05_advanced_analysis_templates.sh 31 trait.txt 8

KMER_LENGTH="$1"
PHENO_FILE="$2"
THREADS="$3"

# 1. Calculate kinship matrix
./bin/emma_kinship_kmers -t kmers_table -k "$KMER_LENGTH" --maf 0.05 > kmers_table.kinship

echo "Kinship matrix created: kmers_table.kinship"

# 2. Convert k-mers table to PLINK binary format (optional)
# ./bin/kmers_table_to_bed -t kmers_table -k "$KMER_LENGTH" -p "$PHENO_FILE" --maf 0.05 --mac 5 -b 10000000 -o output_file
# echo "PLINK binary files created: output_file.bed, .bim, .fam, etc."

# 3. Run k-mers-based GWAS with permutation-based threshold
python2.7 ./kmers_gwas.py --pheno "$PHENO_FILE" --kmers_table kmers_table -l "$KMER_LENGTH" -p "$THREADS" --outdir output_dir --gemma_path ./external_programs/gemma_0_96

echo "GWAS run complete. See output_dir for results."
