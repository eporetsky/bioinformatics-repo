#!/bin/bash
# DADA2 16S Microbiome Analysis Workflow
# This script downloads all samples and runs DADA2 analysis on them together

set -e  # Exit on error

# ============================================================================
# Configuration
# ============================================================================
SAMPLES_TSV="samples.tsv"
FASTQ_DIR="fastq"
DOWNLOAD_DIR="ramdisk"
OUTPUT_DIR="dada2_output"
THREADS=70

# Create directories
mkdir -p "$DOWNLOAD_DIR"
mkdir -p "$FASTQ_DIR"
mkdir -p "$OUTPUT_DIR"

# Step 1: Generate sample list from samples.tsv
if [ ! -f "$SAMPLES_TSV" ]; then
    echo "ERROR: $SAMPLES_TSV not found!" >&2
    exit 1
fi
# Extract the first column minus the header for sample names
SAMPLE_LIST=$(tail -n +2 "$SAMPLES_TSV" | cut -f1)
TOTAL_SAMPLES=$(echo "$SAMPLE_LIST" | wc -l)
echo "Found $TOTAL_SAMPLES samples to process."

CURRENT=0
echo "$SAMPLE_LIST" | while read -r SAMPLE; do
    CURRENT=$((CURRENT + 1))
    echo "[$CURRENT/$TOTAL_SAMPLES] Processing $SAMPLE..."
    FQ1_GZ="$FASTQ_DIR/${SAMPLE}_1.fastq.gz"
    FQ2_GZ="$FASTQ_DIR/${SAMPLE}_2.fastq.gz"
    if [ -f "$FQ1_GZ" ] && [ -f "$FQ2_GZ" ]; then
        echo "  -> Compressed fastq files found in $FASTQ_DIR, skipping download."
        continue
    fi
    # Download if the gzipped fastqs are not present
    echo "  -> Running fasterq-dump into $DOWNLOAD_DIR..."
    fasterq-dump --temp ramdisk/ $SAMPLE --threads $THREADS -O ramdisk/
    
    # Gzip the files
    if [ -f "$DOWNLOAD_DIR/${SAMPLE}_1.fastq" ]; then
        pigz --keep --best "$DOWNLOAD_DIR/${SAMPLE}_1.fastq"
        mv "$DOWNLOAD_DIR/${SAMPLE}_1.fastq.gz" "$FASTQ_DIR/"
    fi
    if [ -f "$DOWNLOAD_DIR/${SAMPLE}_2.fastq" ]; then
        pigz --keep --best "$DOWNLOAD_DIR/${SAMPLE}_2.fastq"
        mv "$DOWNLOAD_DIR/${SAMPLE}_2.fastq.gz" "$FASTQ_DIR/"
    fi
    echo "  -> Download and compression complete for $SAMPLE."
done

echo ""
echo "All samples downloaded successfully!"

# ============================================================================
# Step 2: Quality Control with FastQC (Optional but recommended)
# ============================================================================
echo ""
echo "=== Step 2: Quality Control with FastQC ==="
mkdir -p $OUTPUT_DIR/fastqc_raw

echo "Running FastQC on raw reads..."
fastqc -t $THREADS -o $OUTPUT_DIR/fastqc_raw $FASTQ_DIR/*_1.fastq $FASTQ_DIR/*_2.fastq

echo "Generating MultiQC report..."
multiqc -o $OUTPUT_DIR/fastqc_raw $OUTPUT_DIR/fastqc_raw -n raw_reads_multiqc

echo "FastQC complete. Check $OUTPUT_DIR/fastqc_raw/ for results"

# ============================================================================
# Step 3: Run DADA2 pipeline in R
# ============================================================================
echo ""
echo "=== Step 3: Running DADA2 Pipeline ==="
echo "This step processes all samples together for optimal ASV inference"
echo "This may take a while depending on the number of samples..."
echo ""

Rscript dada2_pipeline.R

# ============================================================================
# Step 4: Clean up (optional)
# ============================================================================
echo ""
echo "=== Step 4: Cleanup ==="
# Skipping cleanup prompt and deletion to preserve all raw and intermediate files for reproducibility and rerunning with different parameters.
echo "Raw and intermediate files have been kept for further analysis."

# ============================================================================
# Done!
# ============================================================================
echo ""
echo "========================================="
echo "DADA2 Analysis Complete!"
echo "========================================="
echo ""
echo "Results are in: $OUTPUT_DIR/"
echo ""
echo "Key output files:"
echo "  - ASVs.fasta: FASTA file of all unique ASV sequences"
echo "  - ASV_counts.tsv: Abundance table (samples x ASVs)"
echo "  - ASV_taxonomy.txt: Taxonomic classification"
echo "  - track_reads.csv: Read counts through each step"
echo "  - phyloseq_object.rds: R phyloseq object for downstream analysis"
echo ""
echo "Next steps:"
echo "  1. Review quality metrics in track_reads.csv"
echo "  2. Load phyloseq_object.rds in R for statistical analysis"
echo "  3. Perform differential abundance testing"
echo "  4. Generate publication-quality figures"
echo ""
