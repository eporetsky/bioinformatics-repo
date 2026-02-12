#!/usr/bin/env bash
# Run Diamond2GO (no InterProScan) on all clean/*.fa files.
# Output: diamond2go/{genome_id}.go.tsv. Skips existing output files.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

D2GO="${SCRIPT_DIR}/Diamond2GO-main/diamond2go"
CLEAN_DIR="${SCRIPT_DIR}/clean"
OUT_DIR="${SCRIPT_DIR}/diamond2go"

if [[ ! -x "$D2GO" && ! -f "$D2GO" ]]; then
  echo "Error: diamond2go not found at $D2GO" >&2
  exit 1
fi

if [[ ! -d "$CLEAN_DIR" ]]; then
  echo "Error: clean directory not found: $CLEAN_DIR" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

for fa in "$CLEAN_DIR"/*.fa; do
  [[ -f "$fa" ]] || continue
  genome_id="${fa##*/}"
  genome_id="${genome_id%.fa}"
  out_tsv="${OUT_DIR}/${genome_id}.go.tsv"

  if [[ -f "$out_tsv" ]]; then
    echo "Skip (exists): $genome_id"
    continue
  fi

  echo "Running: $genome_id"
  "$D2GO" -q "$fa" -b "$out_tsv" -s 12 -r 70 -g 64 -k 32 -m 5
done

echo "Done."
