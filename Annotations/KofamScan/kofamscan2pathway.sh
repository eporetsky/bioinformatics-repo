#!/usr/bin/env bash
# Map KofamScan KO annotations to KEGG pathways
# Downloads KO→pathway and pathway name tables from KEGG REST API,
# then joins them with both clean and clean.hc files to produce gene→pathway TSVs.
#
# Usage: ./kofamscan2pathway.sh [results_dir]
#        Processes all *.clean.tsv and *.clean.hc.tsv in results_dir (default: kofamscan_results/)
#
# Output per input file:
#   {basename}.mapper.txt          - gene_id<TAB>KO from hc (no header, for KEGG Mapper web upload)
#   {basename}.gene2pathway.tsv    - gene→pathway from clean (all predictions)
#   {basename}.gene2pathway.hc.tsv - gene→pathway from clean.hc (high-confidence only)

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${1:-$SCRIPT_DIR/kofamscan_results}"
CACHE_DIR="$SCRIPT_DIR"

KO_PW_FILE="$CACHE_DIR/ko_pathway_link.tsv"
PW_NAMES_FILE="$CACHE_DIR/pathway_names.tsv"

KEGG_API="https://rest.kegg.jp"

if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "Error: results directory not found: $RESULTS_DIR"
    exit 1
fi

# Download KO→pathway link table (cached)
if [[ ! -f "$KO_PW_FILE" ]]; then
    echo "Downloading KO→pathway mapping from KEGG REST API..."
    curl -sf "$KEGG_API/link/pathway/ko" -o "$KO_PW_FILE"
    echo "  $(wc -l < "$KO_PW_FILE") entries"
else
    echo "Using cached $KO_PW_FILE ($(wc -l < "$KO_PW_FILE") entries)"
fi

# Download pathway names (cached)
if [[ ! -f "$PW_NAMES_FILE" ]]; then
    echo "Downloading pathway names from KEGG REST API..."
    curl -sf "$KEGG_API/list/pathway" -o "$PW_NAMES_FILE"
    echo "  $(wc -l < "$PW_NAMES_FILE") pathways"
else
    echo "Using cached $PW_NAMES_FILE ($(wc -l < "$PW_NAMES_FILE") pathways)"
fi

echo ""

shopt -s nullglob
HC_FILES=("$RESULTS_DIR"/*.clean.hc.tsv)
shopt -u nullglob

if [[ ${#HC_FILES[@]} -eq 0 ]]; then
    echo "No *.clean.hc.tsv files found in $RESULTS_DIR"
    echo "Run kofamscan_to_tsv.sh first."
    exit 1
fi

# Maps a two-column (gene TAB KO) file to gene→pathway using the cached KEGG tables
map_to_pathway() {
    local mapper="$1" output="$2"
    awk -F'\t' '
        BEGIN { OFS="\t" }
        FILENAME == ARGV[1] { pathnames[$1] = $2; next }
        FILENAME == ARGV[2] {
            gsub(/^ko:/, "", $1); gsub(/^path:/, "", $2)
            if ($2 ~ /^map/) {
                ko = $1; pw = $2
                ko2pw[ko] = ko2pw[ko] ? ko2pw[ko] SUBSEP pw : pw
            }
            next
        }
        {
            gene = $1; ko = $2
            if (ko in ko2pw) {
                n = split(ko2pw[ko], pws, SUBSEP)
                for (i = 1; i <= n; i++) {
                    pw = pws[i]
                    print gene, ko, pw, (pw in pathnames ? pathnames[pw] : "")
                }
            }
        }
    ' "$PW_NAMES_FILE" "$KO_PW_FILE" "$mapper" > "$output"
}

for hc_file in "${HC_FILES[@]}"; do
    base="${hc_file%.clean.hc.tsv}"
    clean_file="${base}.clean.tsv"
    mapper_hc="${base}.mapper.txt"
    g2p="${base}.gene2pathway.tsv"
    g2p_hc="${base}.gene2pathway.hc.tsv"

    echo "$(basename "$base"):"

    # HC: mapper file (gene_id TAB KO, no header) + gene2pathway
    tail -n +2 "$hc_file" | awk -F'\t' '{print $1 "\t" $2}' > "$mapper_hc"
    map_to_pathway "$mapper_hc" "$g2p_hc"
    echo "  hc:    $(wc -l < "$mapper_hc") genes, $(wc -l < "$g2p_hc") gene-pathway pairs"

    # Clean (all predictions): gene2pathway
    if [[ -f "$clean_file" ]]; then
        mapper_all=$(mktemp)
        tail -n +2 "$clean_file" | awk -F'\t' '{print $2 "\t" $3}' > "$mapper_all"
        map_to_pathway "$mapper_all" "$g2p"
        rm -f "$mapper_all"
        echo "  all:   $(wc -l < "$g2p") gene-pathway pairs"
    else
        echo "  Warning: $clean_file not found, skipping all-predictions mapping"
    fi
done

echo ""
echo "Done."
