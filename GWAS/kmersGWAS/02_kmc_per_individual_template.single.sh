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

# Usage: ./02_kmc_per_individual_template.single.sh <KMER_LENGTH> <THREADS> <sample_name> <SAMPLES_TSV>
# Called by 02_kmc_per_individual_template.sh for each unique sample_name in samples.tsv.

KMER_LENGTH="$1"
THREADS="$2"
sample="$3"
SAMPLES_TSV="${4:-samples.tsv}"

# Resolve samples path if relative
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[[ "$SAMPLES_TSV" != /* ]] && SAMPLES_TSV="$SCRIPT_DIR/$SAMPLES_TSV"

# Create output directory for the sample
outdir="kmc/$sample"
mkdir -p "$outdir"

KMC_BIN="./external_programs/kmc_v3"
STRAND_BIN="./bin/kmers_add_strand_information"

# Build input_files.txt from samples.tsv: for this sample_name, add fq1 and fq2 (if non-empty) for each row.
# Columns: sample_name, library_name, fq1, fq2, sra. Paired: list fq1 then fq2 per row; single: only fq1.
awk -F'\t' -v sample="$sample" '
  NR>1 && $1==sample {
    # fq1 is column 3, fq2 is column 4
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $3);
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $4);
    if ($3 != "") print $3;
    if ($4 != "") print $4
  }
' "$SAMPLES_TSV" > "$outdir/input_files.txt"

if [[ ! -s "$outdir/input_files.txt" ]]; then
  echo "Error: no fastq paths found for sample '$sample' in $SAMPLES_TSV" >&2
  exit 1
fi

# Canonized k-mer counting
"$KMC_BIN" -t"$THREADS" -k"$KMER_LENGTH" -ci2 @"$outdir/input_files.txt" "$outdir/output_kmc_canon" "$outdir" 1> "$outdir/kmc_canon.1" 2> "$outdir/kmc_canon.2"

# Non-canonized k-mer counting
"$KMC_BIN" -t"$THREADS" -k"$KMER_LENGTH" -ci0 -b @"$outdir/input_files.txt" "$outdir/output_kmc_all" "$outdir" 1> "$outdir/kmc_all.1" 2> "$outdir/kmc_all.2"

# Combine canonized and non-canonized
"$STRAND_BIN" -c "$outdir/output_kmc_canon" -n "$outdir/output_kmc_all" -k "$KMER_LENGTH" -o "$outdir/kmers_with_strand"

# Clean up large KMC files
#rm -f "$outdir"/*.kmc*

echo "Done for $sample"
