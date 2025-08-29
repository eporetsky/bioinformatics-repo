#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="S3"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=200GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 2-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

set -e

# Usage: ./03_combine_and_filter_kmers.sh <KMER_LENGTH> <MAC> <PERCENT>
# Example: ./03_combine_and_filter_kmers.sh 31 5 0.2

KMER_LENGTH="$1"
MAC="$2"      # e.g., 5
PERCENT="$3"  # e.g., 0.2

# Generate kmers_list_paths.txt automatically
# kmers_list_paths.txt

# Probably not necessary
rm kmers_list_paths.txt

# Create kmers_list_paths.txt with full paths and sample names
# Each line: full_path_to_kmers_with_strand<tab>sample_name
ls kmc/ | tail -n +2 | awk '{printf "/path/kmer-gwas/kmc/%s/kmers_with_strand\t%s\n", $NF, $NF}' > kmers_list_paths.txt

# Combine and filter
./bin/list_kmers_found_in_multiple_samples -l kmers_list_paths.txt -k "$KMER_LENGTH" --mac "$MAC" -p "$PERCENT" -o kmers_to_use

echo "Filtered k-mers list created: kmers_to_use"
