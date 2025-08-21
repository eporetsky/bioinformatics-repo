#!/bin/bash
set -e

# Usage: ./02_kmc_per_individual_template.sh <KMER_LENGTH> <THREADS>
# Example: ./02_kmc_per_individual_template.sh 31 2

KMER_LENGTH="$1"
THREADS="$2"

FASTQ_DIR="./fastq"
KMC_BIN="./external_programs/kmc_v3"
STRAND_BIN="./bin/kmers_add_strand_information"

# Loop over all unique sample names (assuming files are named SAMPLEID*.fastq)
for fq in "$FASTQ_DIR"/*.fastq; do
    fname=$(basename "$fq")
    sample="${fname%%_*}"
    sample="${sample%%.*}"

    # Create output directory for the sample
    outdir="kmc/$sample"
    mkdir -p "$outdir"

    # Generate input_files.txt for this sample (all matching fastq files in FASTQ_DIR)
    find "$FASTQ_DIR" -type f -name "$sample*.fastq" > "$outdir/input_files.txt"

    # Canonized k-mer counting
    "$KMC_BIN" -t"$THREADS" -k"$KMER_LENGTH" -ci2 @"$outdir/input_files.txt" "$outdir/output_kmc_canon" "$outdir" 1> "$outdir/kmc_canon.1" 2> "$outdir/kmc_canon.2"

    # Non-canonized k-mer counting
    "$KMC_BIN" -t"$THREADS" -k"$KMER_LENGTH" -ci0 -b @"$outdir/input_files.txt" "$outdir/output_kmc_all" "$outdir" 1> "$outdir/kmc_all.1" 2> "$outdir/kmc_all.2"

    # Combine canonized and non-canonized
    "$STRAND_BIN" -c "$outdir/output_kmc_canon" -n "$outdir/output_kmc_all" -k "$KMER_LENGTH" -o "$outdir/kmers_with_strand"

    # Clean up large KMC files
    rm -f "$outdir"/*.kmc*

    echo "Done for $sample"
done
