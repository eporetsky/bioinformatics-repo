#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="S2"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=200GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 2-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

date

set -e

# Usage: ./02_kmc_per_individual_template.sh <KMER_LENGTH> <THREADS> [SAMPLES_TSV] [bash|sbatch]
# Example: ./02_kmc_per_individual_template.sh 31 2
# Example: ./02_kmc_per_individual_template.sh 31 2 samples.tsv sbatch
# Run mode: "bash" runs each sample in the current shell; "sbatch" submits each sample as a batch job (default).

KMER_LENGTH="$1"
THREADS="$2"
SAMPLES_TSV="${3:-samples.tsv}"
RUN_MODE="${4:-sbatch}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SINGLE_SCRIPT="$SCRIPT_DIR/02_kmc_per_individual_template.single.sh"

# Resolve samples path if relative
[[ "$SAMPLES_TSV" != /* ]] && SAMPLES_TSV="$SCRIPT_DIR/$SAMPLES_TSV"

if [[ ! -f "$SAMPLES_TSV" ]]; then
  echo "Error: samples TSV not found: $SAMPLES_TSV" >&2
  exit 1
fi

if [[ "$RUN_MODE" != "bash" && "$RUN_MODE" != "sbatch" ]]; then
  echo "Error: RUN_MODE must be 'bash' or 'sbatch', got: $RUN_MODE" >&2
  exit 1
fi

# Get unique sample_name values (column 1), skip header.
# In bash mode we don't use set -e for the per-sample run so one failure doesn't stop the rest.
FAILED=""
RUN_LOG=""
while IFS= read -r sample; do
  if [[ -z "$sample" ]]; then continue; fi
  if [[ "$RUN_MODE" == "sbatch" ]]; then
    sbatch "$SINGLE_SCRIPT" "$KMER_LENGTH" "$THREADS" "$sample" "$SAMPLES_TSV"
  else
    run_out="$(mktemp)"
    if ! bash "$SINGLE_SCRIPT" "$KMER_LENGTH" "$THREADS" "$sample" "$SAMPLES_TSV" > "$run_out" 2>&1; then
      FAILED="${FAILED:+$FAILED }$sample"
      RUN_LOG="${RUN_LOG}$(printf '\n--- %s failed ---\n%s\n' "$sample" "$(cat "$run_out")")"
    fi
    rm -f "$run_out"
  fi
done < <(awk -F'\t' 'NR>1 {print $1}' "$SAMPLES_TSV" | sort -u)

if [[ -n "$FAILED" ]]; then
  echo "Warning: the following samples failed: $FAILED" >&2
  echo "$RUN_LOG" >&2
  exit 1
fi
