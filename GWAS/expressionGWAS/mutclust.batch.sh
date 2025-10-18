#!/bin/bash

# Batch script to run MutClust MR and CLS on all CPM files
# Usage: ./mutclust.batch.sh

set -e  # Exit on error

# Find all .cpm.tsv.gz files in data/*/cpm/ directories
for cpm_file in data/*.tsv; do
    # Check if file exists (in case glob doesn't match anything)
    if [ ! -f "$cpm_file" ]; then
        echo "No CPM files found matching pattern: data/*tsv"
        exit 1
    fi
    
    # Extract genotype and accession from path
    # Example: data/AtCol0/cpm/PRJEB9918.cpm.tsv.gz
    condition=$(basename "$cpm_file" .tsv)
    
    echo "=========================================="
    echo "Processing: $condition"
    echo "=========================================="
    
    # Define output paths
    mr_output="MR/${condition}.mr.tsv.gz"
    cls_output="CLS/${condition}.cls5.tsv"
    
    # Create output directories if they don't exist
    mkdir -p "MR"
    mkdir -p "CLS"
    
    # Check the number of columns in the RNAseq input (excluding first field)
    n_cols=$(zcat "$cpm_file" | head -n1 | awk -F'\t' '{print NF-1}')
    if [ "$n_cols" -lt 12 ]; then
        echo "  -> Skipping: Input $cpm_file has only $n_cols samples (<12)."
        # Remove outputs if they exist
        if [ -f "$mr_output" ]; then
            echo "    -> Deleting $mr_output due to low column count."
            rm -f "$mr_output"
        fi
        if [ -f "$cls_output" ]; then
            echo "    -> Deleting $cls_output due to low column count."
            rm -f "$cls_output"
        fi
        continue
    fi

    # Step 1: Calculate mutual rank (Skip if already exists)
    echo "[1/2] Calculating mutual rank..."
    if [ -f "$mr_output" ]; then
        echo "  -> Skipping: $mr_output already exists"
    else
        mutclust mr --input "$cpm_file" --output "$mr_output" -m 200 -t 70
        echo "  -> Saved: $mr_output"
    fi
    
    # Step 2: Run clustering
    echo "[2/2] Running clustering..."
    if [ -f "$cls_output" ]; then
        echo "  -> Skipping: $cls_output already exists"
    else
        mutclust cls --input "$mr_output" --output "$cls_output" -e 5
        echo "  -> Saved: $cls_output"
    fi
    
    echo "Completed: $condition"
    echo ""
done

echo "=========================================="
echo "All files processed successfully!"
echo "=========================================="

