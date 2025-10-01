#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="S2"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=200GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 2-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

date

set -e

# Usage: ./02_kmc_per_individual_template.sh <KMER_LENGTH> <THREADS>
# Example: ./02_kmc_per_individual_template.sh 31 2

KMER_LENGTH="$1"
THREADS="$2"

FASTQ_DIR="/path/kmer-gwas/fastq"
KMC_BIN="./external_programs/kmc_v3"
STRAND_BIN="./bin/kmers_add_strand_information"

# Loop over all unique sample names (assuming files are named SAMPLEID*.fastq)
for fq in "$FASTQ_DIR"/*.fastq; do
    fname=$(basename "$fq")
    sample="${fname%%_*}"
    sample="${sample%%.*}"

    sbatch 02_kmc_per_individual_template.single.sh "$KMER_LENGTH" "$THREADS" "$sample"
done