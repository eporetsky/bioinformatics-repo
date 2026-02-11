# **WIP:** kmersGWAS Pipeline Setup and Usage

This repository contains scripts and environment setup to help with running the k-mers-based GWAS analysis using the [kmersGWAS package](https://github.com/voichek/kmersGWAS/tree/master/src).

## 1. Environment Setup

1. **Create the conda environment:**
   ```bash
   conda env create -f environment.yml
   conda activate kmersgwas
   ```
2. **Install the R package `mvnpermute` (not available via conda):**
   ```bash
   R
   > install.packages("mvnpermute")
   # Then quit R
   ```

## 2. Download and Set Up kmersGWAS

Run the setup script to download kmersGWAS, KMC, and GEMMA:
```bash
bash 01_download_kmersGWAS.sh
```
This will download all the repository into the current directory with all required binaries.

## 3. Prepare Your Data
- **Sample manifest:** Create a tab-separated file `samples.tsv` with columns: `sample_name`, `library_name`, `fq1`, `fq2`, `sra`. For each row, set `fq1` to the path of the first (or only) FASTQ file and `fq2` to the second file for paired-end data, or leave `fq2` empty for single-end. One `sample_name` can appear on multiple rows (e.g. multiple runs or libraries); all listed FASTQs for that sample will be combined into one k-mer count.
- Ensure your phenotype table is available (e.g., `trait.txt`).
- **Test data:** The example test files (e.g. `test_reads/`, `samples.tsv`, `disease.pheno`) were downloaded from [kGWASflow](https://github.com/akcorut/kGWASflow/).

## 4. Run the Pipeline

### Step 2: KMC Counting for Each Individual
This step reads `samples.tsv`, collects all `fq1` (and `fq2` when present) paths per unique `sample_name`, and runs KMC for each sample. Per-sample input is written to `kmc/<sample_name>/input_files.txt` in the format expected by KMC (one path per line, R1/R2 pairs kept in order).

**Usage:**
```bash
bash 02_kmc_per_individual_template.sh <KMER_LENGTH> <THREADS> [SAMPLES_TSV] [bash|sbatch]
```

**Examples:**
```bash
# Submit one batch job per sample (default)
bash 02_kmc_per_individual_template.sh 31 2

# Use a custom samples file and run mode
bash 02_kmc_per_individual_template.sh 31 2 samples.tsv sbatch

# Run all samples in the current shell (no job submission)
bash 02_kmc_per_individual_template.sh 31 2 samples.tsv bash
```

- **Arguments:** `KMER_LENGTH` (e.g. 31), `THREADS` (e.g. 2), optional `SAMPLES_TSV` (default: `samples.tsv`), optional `bash` or `sbatch` (default: `sbatch`).
- **Output:** For each sample, a directory `kmc/<sample_name>/` containing `kmers_with_strand`, `input_files.txt`, and logs.
- **If you see fewer samples in `kmc/` than in your TSV:** In `bash` mode the launcher now runs all samples even if one fails; it reports which samples failed at the end. Fix any missing FASTQ paths or errors for those samples and re-run (you can re-run only the failed samples if needed).

### Step 3: Combine and Filter k-mers
The script auto-generates `kmers_list_paths.txt` from existing `kmc/<sample_name>/kmers_with_strand` files (from step 02), then combines and filters k-mers:
```bash
bash 03_combine_and_filter_kmers.sh 31 5 0.2
```
- Arguments: `<KMER_LENGTH> <MAC> <PERCENT>`
- Requires step 02 to be complete (all samples in `kmc/`).

### Step 4: Build the k-mers Table
```bash
bash 04_build_kmers_table.sh 31
```
- Arguments: `<KMER_LENGTH>`

### Step 5: Build Kinship Matrix
```bash
bash 05_build_kinship_matrix.sh 31
```
- Arguments: `<KMER_LENGTH>`
- Requires step 04 (kmers_table).

### Step 6: Run GWAS
```bash
bash 06_run_gwas.sh 31 disease.pheno 8
```
- Arguments: `<KMER_LENGTH> [PHENO_FILE] [THREADS]` — default phenotype: `disease.pheno`, default threads: 8.
- Requires step 05 (kmers_table.kinship).

## Notes
- Adjust script parameters as needed for your dataset.
- Phenotype file is tab-separated with “accession_id” and “phenotype_value” headers
- Cite the original work: Voichek, Y., & Weigel, D. (2020). Identifying genetic variants underlying phenotypic variation in plants without complete genomes. Nature genetics, 52(5), 534-540.

## Troubleshooting

**`libgfortran.so.3: cannot open shared object file` (Step 6 / GEMMA)**  
The bundled GEMMA binary needs `libgfortran.so.3`. We keep this out of `environment.yml` to avoid dependency conflicts with R and other packages. Install it only when you run step 6:

- **Conda (e.g. on Ubuntu 22.04+ where system libgfortran3 is not available):**  
  `conda install -c conda-forge libgfortran=3.0.0`  
  Then run step 6 with the env activated; the script sets `LD_LIBRARY_PATH` to the env’s `lib`.
- **System (older Ubuntu/Debian):**  
  `sudo apt install libgfortran3`  
  if the package is available.

**`SafetyError: The package for r-base ... appears to be corrupted`**  
A cached R (or r-base) package in conda’s package cache is corrupted. Remove that package from the cache so conda can re-download it. The error message shows the path, e.g. `/path/to/miniconda3/pkgs/r-base-4.0.3-ha43b4e8_3`. Run:

```bash
rm -rf /path/to/miniconda3/pkgs/r-base-4.0.3-ha43b4e8_3
```

Use the exact path from your error (replace `/path/to/miniconda3` with your conda install path, e.g. `$HOME/miniconda3`). Then retry:

```bash
conda env create -f environment.yml
```

If another package is reported as corrupted, remove that package’s folder under `pkgs/` in the same way and retry.

**`No such file or directory: 'output_dir/pheno.phenotypes_and_permutations'`**  
The R step (phenotype transformation) did not produce this file. Check `output_dir/phenotypes_transformation_permutation.log`. If it says *"Kinship matrix is not positive semi-definite"*, the kinship matrix from step 05 is not valid for the R code; try a higher MAF in step 05 (e.g. `--maf 0.1`) or use more samples. If the log shows another R error, fix that (e.g. missing R packages).

**Step 6 or R stopped working after changing the environment**  
If you added `libgfortran` (or other packages) to `environment.yml` and the pipeline or R started failing, try recreating the env *without* libgfortran: remove it from `environment.yml`, then `conda env remove -n kmersgwas`, `conda env create -f environment.yml`, reinstall R’s `mvnpermute`, and run the pipeline again. Install libgfortran only when you need to run step 6 (see above).