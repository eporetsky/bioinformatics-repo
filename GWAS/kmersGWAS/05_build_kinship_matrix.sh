#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="S5"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=200GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 2-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

set -e

# Usage: ./05_build_kinship_matrix.sh <KMER_LENGTH>
# Example: ./05_build_kinship_matrix.sh 31

KMER_LENGTH="$1"

# Calculate kinship matrix
./bin/emma_kinship_kmers -t kmers_table -k "$KMER_LENGTH" --maf 0.05 > kmers_table.kinship

echo "Kinship matrix created: kmers_table.kinship"