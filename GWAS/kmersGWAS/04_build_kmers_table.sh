#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="S4"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=200GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 2-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

set -e

# Usage: ./04_build_kmers_table.sh <KMER_LENGTH>
# Example: ./04_build_kmers_table.sh 31

KMER_LENGTH="$1"

# Requires kmers_list_paths.txt and kmers_to_use from step 03
if [[ ! -f kmers_list_paths.txt ]]; then
  echo "Error: kmers_list_paths.txt not found. Run step 03 first." >&2
  exit 1
fi
if [[ ! -f kmers_to_use ]]; then
  echo "Error: kmers_to_use not found. Run step 03 first." >&2
  exit 1
fi

./bin/build_kmers_table -l kmers_list_paths.txt -k "$KMER_LENGTH" -a kmers_to_use -o kmers_table

echo "k-mers table created: kmers_table.table and kmers_table.names"
