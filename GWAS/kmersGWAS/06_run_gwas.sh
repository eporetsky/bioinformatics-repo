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

# I had an unresolved issue getting this to work on the cluster (error message below), but it seems to work fine on workstation
# ./external_programs/gemma_0_96: error while loading shared libraries: libgfortran.so.3: cannot open shared object file: No such file or directory

# Usage: ./06_run_gwas.sh <KMER_LENGTH> <PHENO_FILE> <THREADS>
# Example: ./06_run_gwas.sh 31 trait.txt 70

KMER_LENGTH="$1"
PHENO_FILE="$2"
THREADS="$3"

# 3. Run k-mers-based GWAS with permutation-based threshold
rm -r output_dir
python2.7 ./kmers_gwas.py --pheno "$PHENO_FILE" --kmers_table kmers_table -l "$KMER_LENGTH" -p "$THREADS" --outdir output_dir --gemma_path ./external_programs/gemma_0_96

echo "GWAS run complete. See output_dir for results."
