#!/bin/bash
#SBATCH --account=account_name
#SBATCH --partition=partition_name
#SBATCH --job-name="S5"
#SBATCH -N1
#SBATCH -n1
#SBATCH --mem=200GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH -t 2-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

set -e

# Usage: bash 06_run_gwas.sh <KMER_LENGTH> <PHENO_FILE> <THREADS>
# Example: bash 06_run_gwas.sh 31 disease.pheno 8
#
# If you see "libgfortran.so.3: cannot open shared object file", install it:
#   conda install -c conda-forge libgfortran=3.0.0
# (On older Ubuntu/Debian you can use: sudo apt install libgfortran3)

KMER_LENGTH="$1"
PHENO_FILE="${2:-disease.pheno}"   # default disease.pheno for testing
THREADS="${3:-8}"

if [[ ! -f "$PHENO_FILE" ]]; then
  echo "Error: phenotype file not found: $PHENO_FILE" >&2
  exit 1
fi

# GEMMA needs libgfortran.so.3 (not .so.5). Find it and set LD_LIBRARY_PATH.
GEMMA_LIB=""
for dir in "${CONDA_PREFIX:-}/lib" /usr/lib/x86_64-linux-gnu /usr/lib external_programs/lib; do
  [[ -z "$dir" || ! -d "$dir" ]] && continue
  if [[ -f "$dir/libgfortran.so.3" || -f "$dir/libgfortran.so.3.0.0" ]]; then
    GEMMA_LIB="$dir"
    break
  fi
done
if [[ -n "$GEMMA_LIB" ]]; then
  export LD_LIBRARY_PATH="${GEMMA_LIB}${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
else
  echo "Error: libgfortran.so.3 not found (required by GEMMA). Install with:" >&2
  echo "  conda install -c conda-forge libgfortran=3.0.0" >&2
  echo "Then run this script again with your conda env activated." >&2
  exit 1
fi

# Run k-mers-based GWAS with permutation-based threshold
rm -rf output_dir
EXTRA_ARGS=()
[[ "$4" == "--dont_remove_intermediates" ]] && EXTRA_ARGS+=(--dont_remove_intermediates)
python2.7 ./kmers_gwas.py --pheno "$PHENO_FILE" --kmers_table kmers_table -l "$KMER_LENGTH" -p "$THREADS" --outdir output_dir --gemma_path ./external_programs/gemma_0_96 "${EXTRA_ARGS[@]}"

echo "GWAS run complete. See output_dir for results."
