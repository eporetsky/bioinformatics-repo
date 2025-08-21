#!/bin/bash
set -e

# Usage: ./04_build_kmers_table.sh <KMER_LENGTH>
# Example: ./04_build_kmers_table.sh 31

KMER_LENGTH="$1"

# Requires kmers_list_paths.txt and kmers_to_use from previous step

./bin/build_kmers_table -l kmers_list_paths.txt -k "$KMER_LENGTH" -a kmers_to_use -o kmers_table

echo "k-mers table created: kmers_table.table and kmers_table.names"
