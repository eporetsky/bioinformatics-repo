#!/usr/bin/env bash
# Run KofamScan on all FASTA files in the fasta directory
# Requires: conda env 'kegg' activated, setup_kofamscan.sh completed
#
# NOTE: KofamScan expects PROTEIN sequences (amino acids), not nucleotides.
# Nucleotide FASTA files will produce errors or incorrect results.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTA_DIR="${1:-$SCRIPT_DIR/fasta}"
OUTPUT_DIR="${2:-$SCRIPT_DIR/kofamscan_results}"
CPU="${3:-$(nproc)}"
KOFAM_BIN="$SCRIPT_DIR/kofam_scan-1.3.0"
CONFIG="$KOFAM_BIN/config.yml"

# Check prerequisites
if [[ ! -f "$KOFAM_BIN/exec_annotation" ]]; then
    echo "Error: exec_annotation not found at $KOFAM_BIN"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "Error: config.yml not found. Run setup_kofamscan.sh first:"
    echo "  conda activate kegg"
    echo "  ./setup_kofamscan.sh"
    exit 1
fi

if [[ ! -d "$FASTA_DIR" ]]; then
    echo "Error: FASTA directory not found: $FASTA_DIR"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Find all FASTA files (.fa, .fasta, .faa, .fna)
shopt -s nullglob
FASTA_FILES=("$FASTA_DIR"/*.fa "$FASTA_DIR"/*.fasta "$FASTA_DIR"/*.faa "$FASTA_DIR"/*.fna)
shopt -u nullglob

# Deduplicate (in case of overlapping globs)
readarray -t FASTA_FILES < <(printf '%s\n' "${FASTA_FILES[@]}" | sort -u)

if [[ ${#FASTA_FILES[@]} -eq 0 ]]; then
    echo "No FASTA files found in $FASTA_DIR"
    exit 0
fi

echo "Found ${#FASTA_FILES[@]} FASTA file(s)"
echo "Output directory: $OUTPUT_DIR"
echo "Using $CPU CPU cores"
echo ""

for fasta in "${FASTA_FILES[@]}"; do
    basename=$(basename "$fasta")
    name="${basename%.*}"
    output_tsv="$OUTPUT_DIR/${name}_kofamscan.tsv"

    echo "Processing: $basename"
    if "$KOFAM_BIN/exec_annotation" \
        -c "$CONFIG" \
        --cpu "$CPU" \
        -o "$output_tsv" \
        -f detail-tsv \
        -E 0.00001 \
        --report-unannotated \
        "$fasta"; then
        echo "  -> $output_tsv"
    else
        echo "  Warning: KofamScan failed for $basename (KofamScan requires protein sequences, not nucleotide)"
    fi
done

echo ""
echo "Done. Results in $OUTPUT_DIR"
