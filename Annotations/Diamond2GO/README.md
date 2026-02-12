# Diamond2GO batch pipeline

Runs [Diamond2GO](https://github.com/rhysf/Diamond2GO) on many FASTA files to assign Gene Ontology (GO) terms. Uses DIAMOND only (no InterProScan).

## Setup

1. **Conda environment**

   ```bash
   conda env create -f environment.yml
   conda activate diamond2go
   ```

2. **Reference database**

   Put the Diamond2GO DIAMOND database in `Diamond2GO-main/resources/` (e.g. `nr_clean_d2go_20250812_c95.faa.dmnd`). It will download automatically the first time you run the script. To skip the MD5 checks and save time, add a `return 1;` at the start of the `check_or_download_db` function in `Diamond2GO-main/util/Diamond2go.pl`.

## Input / output

- **Input:** one FASTA per genome in `clean/`, named `{genome_id}.fa` (e.g. `clean/AtCol-0.fa`).
- **Output:** one TSV per genome in `diamond2go/`, named `{genome_id}.go.tsv`, with columns: gene_id, gene_name, species, e-value, GO-term, evidence_code, qualifier, category.

## Run

From the repo root:

```bash
conda activate diamond2go
./run_diamond2go_batch.sh
```

The script loops over every `clean/*.fa`, runs Diamond2GO (steps 1 and 2 only), and writes `diamond2go/{genome_id}.go.tsv`. Existing output files are skipped.

## Options

Pipeline options (e.g. `-r`, `-g`, `-k`, `-m`) are set in `run_diamond2go_batch.sh`. Edit that script to change threads, block size, or max target sequences.
