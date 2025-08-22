#!/bin/bash
set -e

# Usage: ./03_combine_and_filter_kmers.sh <KMER_LENGTH> <MAC> <PERCENT>
# Example: ./03_combine_and_filter_kmers.sh 31 5 0.2

KMER_LENGTH="$1"
MAC="$2"      # e.g., 5
PERCENT="$3"  # e.g., 0.2


# Create kmers_list_paths.txt with full paths and sample names
# Each line: full_path_to_kmers_with_strand<tab>sample_name
ls kmc/ | tail -n +2 | awk '{printf "kmc/%s/kmers_with_strand\t%s\n", $NF, $NF}' > kmers_list_paths.txt

# Combine and filter
./bin/list_kmers_found_in_multiple_samples -l kmers_list_paths.txt -k "$KMER_LENGTH" --mac "$MAC" -p "$PERCENT" -o kmers_to_use

echo "Filtered k-mers list created: kmers_to_use"
