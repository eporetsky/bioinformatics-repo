#!/bin/bash
set -e

# Usage: ./03_combine_and_filter_kmers.sh <KMER_LENGTH> <MAC> <PERCENT>
# Example: ./03_combine_and_filter_kmers.sh 31 5 0.2

KMER_LENGTH="$1"
MAC="$2"      # e.g., 5
PERCENT="$3"  # e.g., 0.2

# Generate kmers_list_paths.txt automatically
> kmers_list_paths.txt
for d in kmc/*/; do
    sample=$(basename "$d")
    if [ -f "$d/kmers_with_strand" ]; then
        echo "$d/kmers_with_strand\t$sample" >> kmers_list_paths.txt
    fi
done

# Combine and filter
./bin/list_kmers_found_in_multiple_samples -l kmers_list_paths.txt -k "$KMER_LENGTH" --mac "$MAC" -p "$PERCENT" -o kmers_to_use

echo "Filtered k-mers list created: kmers_to_use"
