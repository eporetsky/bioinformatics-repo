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

# Usage: bash 03_combine_and_filter_kmers.sh <KMER_LENGTH> <MAC> <PERCENT>
# Example: bash 03_combine_and_filter_kmers.sh 31 5 0.2

KMER_LENGTH="$1"
MAC="$2"      # e.g., 5
PERCENT="$3"  # e.g., 0.2

if [[ -z "$KMER_LENGTH" || -z "$MAC" || -z "$PERCENT" ]]; then
  echo "Usage: bash 03_combine_and_filter_kmers.sh <KMER_LENGTH> <MAC> <PERCENT>" >&2
  echo "Example: bash 03_combine_and_filter_kmers.sh 31 5 0.2" >&2
  exit 1
fi

# Generate kmers_list_paths.txt from kmc/<sample_name>/kmers_with_strand (output of step 02)
# Each line: full_path_to_kmers_with_strand<tab>sample_name
rm -f kmers_list_paths.txt
base="$(pwd)"
for dir in kmc/*/; do
  [[ -d "$dir" ]] || continue
  sample="$(basename "$dir")"
  f="kmc/$sample/kmers_with_strand"
  if [[ -f "$f" ]]; then
    printf "%s/%s\t%s\n" "$base" "$f" "$sample" >> kmers_list_paths.txt
  fi
done
if [[ ! -s kmers_list_paths.txt ]]; then
  echo "Error: no kmc/<sample>/kmers_with_strand files found. Run step 02 first." >&2
  exit 1
fi

# Combine and filter
./bin/list_kmers_found_in_multiple_samples -l kmers_list_paths.txt -k "$KMER_LENGTH" --mac "$MAC" -p "$PERCENT" -o kmers_to_use

echo "Filtered k-mers list created: kmers_to_use"
