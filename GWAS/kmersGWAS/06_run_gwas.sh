#!/bin/bash
set -e

# Usage: ./06_run_gwas.sh <KMER_LENGTH> <PHENO_FILE> <THREADS>
# Example: ./06_run_gwas.sh 31 trait.txt 70

KMER_LENGTH="$1"
PHENO_FILE="$2"
THREADS="$3"

# 3. Run k-mers-based GWAS with permutation-based threshold

rm -r output_dir
python2.7 ./kmers_gwas.py --pheno "$PHENO_FILE" --kmers_table kmers_table -l "$KMER_LENGTH" -p "$THREADS" --outdir output_dir --gemma_path ./external_programs/gemma_0_96

echo "GWAS run complete. See output_dir for results."