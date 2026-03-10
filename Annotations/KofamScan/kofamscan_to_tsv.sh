#!/usr/bin/env bash
# Convert KofamScan detail-tsv output to proper TSV with clean headers
# Replaces #-commented headers with proper column names; converts * to hc=y
#
# Usage: ./kofamscan_to_tsv.sh [dir|input.tsv] [output.tsv]
#        Default: ./kofamscan_to_tsv.sh kofamscan_results/
#        Converts all *_kofamscan.tsv in directory to *.clean.tsv and *.clean.hc.tsv

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT="${1:-$SCRIPT_DIR/kofamscan_results}"
OUTPUT="${2:-}"

# Headers
HEADER="hc	gene_name	KO	threshold	score	e_value	ko_definition"
HEADER_HC="gene_name	KO	threshold	score	e_value	ko_definition"

convert() {
    local file="${1:-$INPUT}"
    echo "$HEADER"
    awk -F'\t' '
        /^#/ { next }
        NF == 0 { next }
        {
            # Col1: * or empty, Col2: gene, Col3: KO, Col4: thrshld, Col5: score, Col6: E-value, Col7: definition
            col1 = ($1 == "*") ? "y" : "n"
            gene = $2
            ko = $3
            thr = $4
            score = $5
            eval = $6
            def = $7
            # Strip quotes from definition
            gsub(/^"|"$/, "", def)
            print col1 "\t" gene "\t" ko "\t" thr "\t" score "\t" eval "\t" def
        }
    ' "$file"
}

# High-confidence only: filters to hc=y, omits hc column
convert_hc() {
    local file="${1:-$INPUT}"
    echo "$HEADER_HC"
    awk -F'\t' '
        /^#/ { next }
        NF == 0 { next }
        $1 != "*" { next }
        {
            gene = $2; ko = $3; thr = $4; score = $5; eval = $6; def = $7
            gsub(/^"|"$/, "", def)
            print gene "\t" ko "\t" thr "\t" score "\t" eval "\t" def
        }
    ' "$file"
}

# Batch mode: convert all *_kofamscan.tsv in directory
if [[ -d "$INPUT" ]]; then
    for f in "$INPUT"/*_kofamscan.tsv; do
        [[ -f "$f" ]] || continue
        base="${f%_kofamscan.tsv}"
        convert "$f" > "${base}_kofamscan.clean.tsv"
        convert_hc "$f" > "${base}_kofamscan.clean.hc.tsv"
        echo "Wrote ${base}_kofamscan.clean.tsv"
        echo "Wrote ${base}_kofamscan.clean.hc.tsv"
    done
    exit 0
fi

if [[ -n "$OUTPUT" ]]; then
    convert > "$OUTPUT"
    convert_hc "$INPUT" > "${OUTPUT%.tsv}.hc.tsv"
    echo "Wrote $OUTPUT"
    echo "Wrote ${OUTPUT%.tsv}.hc.tsv"
else
    convert
fi
