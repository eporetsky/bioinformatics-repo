#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Build a combined, non-redundant FASTA from many {accession}/dada2_output/ASVs.fasta
and create a vsearch UDB database for closed-reference mapping.

Usage:
  ./vsearch_build_closed_reference_db.sh \
    --root /path/to/parent/containing/accessions \
    --out  /path/to/output_dir \
    --id   1.0 \
    --threads 8

Defaults:
  --root    : current working directory
  --out     : <root>/closed_reference_db
  --id      : 1.0   (set 0.99 for 99% identity clustering, etc.)
  --threads : (auto) nproc

Outputs (in --out):
  - ASVs_all.fasta                 (all ASVs, headers prefixed with accession)
  - ASVs_all.derep.fasta           (dereplicated)
  - ASVs_closedref.fasta           (centroids after clustering at --id)
  - ASVs_closedref.udb             (vsearch database)
  - derep.uc, cluster.uc           (membership mappings)
  - inputs.tsv                     (which ASVs.fasta were used)

Notes:
  - This script looks for files matching */dada2_output/ASVs.fasta under --root.
  - Accessions are inferred as the directory directly containing dada2_output/.
EOF
}

ROOT="$(pwd)"
OUT=""
ID="0.97"
THREADS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --id) ID="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "${OUT}" ]]; then
  OUT="${ROOT%/}/closed_reference_db"
fi

if [[ -z "${THREADS}" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    THREADS="$(nproc)"
  else
    THREADS="1"
  fi
fi

if ! command -v vsearch >/dev/null 2>&1; then
  echo "ERROR: vsearch not found in PATH. Install it first (e.g. conda install -c bioconda vsearch)." >&2
  exit 1
fi

mkdir -p "${OUT}"

INPUTS_TSV="${OUT%/}/inputs.tsv"
ALL_FASTA="${OUT%/}/ASVs_all.fasta"
DEREP_FASTA="${OUT%/}/ASVs_all.derep.fasta"
CLOSEDREF_FASTA="${OUT%/}/ASVs_closedref.fasta"
UDB="${OUT%/}/ASVs_closedref.udb"
DEREP_UC="${OUT%/}/derep.uc"
CLUSTER_UC="${OUT%/}/cluster.uc"

rm -f "${ALL_FASTA}" "${INPUTS_TSV}"
printf "accession\tasvs_fasta\n" > "${INPUTS_TSV}"

echo "Parsing ${ROOT}/*/dada2_output/ASVs.fasta (one per accession)"

mapfile -t FASTAS < <(
  for d in "${ROOT}"/*; do
    if [[ -d "$d" && -e "$d/dada2_output/ASVs.fasta" ]]; then
      echo "$d/dada2_output/ASVs.fasta"
    fi
  done
)

if [[ "${#FASTAS[@]}" -eq 0 ]]; then
  echo "ERROR: Found 0 ASVs.fasta files under top-level dirs in ${ROOT}" >&2
  exit 1
fi

for f in "${FASTAS[@]}"; do
  # accession is the directory that contains dada2_output/
  # .../<accession>/dada2_output/ASVs.fasta
  accession="$(basename "$(dirname "$(dirname "$f")")")"

  printf "%s\t%s\n" "${accession}" "${f}" >> "${INPUTS_TSV}"

  # Prefix FASTA headers with accession to avoid collisions across projects.
  # Keep only the first token of the original header (strip after whitespace).
  awk -v acc="${accession}" '
    BEGIN { OFS="" }
    /^>/ {
      h = $0
      sub(/^>/, "", h)
      split(h, a, /[ \t]/)
      id = a[1]
      print ">", acc, "|", id
      next
    }
    { print $0 }
  ' "${f}" >> "${ALL_FASTA}"
done

echo "Wrote combined FASTA: ${ALL_FASTA}"

echo "Dereplicating (full length)..."
vsearch \
  --derep_fulllength "${ALL_FASTA}" \
  --output "${DEREP_FASTA}" \
  --uc "${DEREP_UC}" \
  --sizeout \
  --fasta_width 0 \
  --threads "${THREADS}"

echo "Clustering at identity ${ID}..."
vsearch \
  --cluster_fast "${DEREP_FASTA}" \
  --id "${ID}" \
  --centroids "${CLOSEDREF_FASTA}" \
  --uc "${CLUSTER_UC}" \
  --sizeout \
  --fasta_width 0 \
  --threads "${THREADS}"

echo "Building vsearch database (UDB)..."
vsearch \
  --makeudb_usearch "${CLOSEDREF_FASTA}" \
  --output "${UDB}"

cat <<EOF
Done.

Inputs:              ${INPUTS_TSV}
All ASVs:            ${ALL_FASTA}
Dereplicated:        ${DEREP_FASTA}
Closed-ref FASTA:    ${CLOSEDREF_FASTA}
Closed-ref DB (UDB): ${UDB}
Mappings:            ${DEREP_UC}, ${CLUSTER_UC}
EOF





