# KofamScan Analysis Workflow

Run KofamScan (KEGG Orthology annotation) on all FASTA files in the `fasta/` directory, convert results to clean TSVs, and map KO IDs to KEGG pathways.

## Quick Start

```bash
# 1. Create and activate the conda environment
conda env create -f environment.yml
conda activate kofamscan

# 2. One-time setup: downloads kofam_scan tool, KOfam db (~4GB) -> db/, kofam_scan-1.3.0/
./setup_kofamscan.sh

# 3. Run KofamScan on fasta/*.fa -> kofamscan_results/*_kofamscan.tsv
./run_kofamscan.sh

# 4. Convert raw tsv -> kofamscan_results/*.clean.tsv + *.clean.hc.tsv
./kofamscan_to_tsv.sh

# 5. Map KO to pathways -> kofamscan_results/*.gene2pathway.tsv + *.gene2pathway.hc.tsv
./kofamscan2pathway.sh
```

## Requirements

- **Conda** (Anaconda or Miniconda)
- **Protein sequences only**: KofamScan expects amino acid FASTA files, not nucleotide sequences. Files like `*.protein.fa` are appropriate; nucleotide assemblies will fail.

## Scripts

| Script | Purpose |
|--------|---------|
| `environment.yml` | Conda env definition (ruby, hmmer, parallel) |
| `setup_kofamscan.sh` | Downloads KOfam database, creates `config.yml` |
| `run_kofamscan.sh` | Runs KofamScan on all FASTA files in `fasta/` |
| `kofamscan_to_tsv.sh` | Converts raw KofamScan output to clean TSVs |
| `kofamscan2pathway.sh` | Maps KO IDs to KEGG pathways via REST API |

## Step-by-step

### 1. Setup

```bash
conda env create -f environment.yml
conda activate kegg
./setup_kofamscan.sh
```

This downloads the KOfam HMM profiles and `ko_list` database (~4 GB), and writes a `config.yml` pointing to your conda-installed `hmmsearch` and `parallel`.

### 2. Run KofamScan

```bash
./run_kofamscan.sh
```

Scans all `.fa`, `.fasta`, `.faa`, `.fna` files in `fasta/`. Results go to `kofamscan_results/`.

**Parallelization:** Uses all CPU cores by default. Override with the third argument:

```bash
./run_kofamscan.sh fasta/ kofamscan_results/ 16   # use 16 cores
```

### 3. Convert to clean TSV

```bash
./kofamscan_to_tsv.sh
```

Defaults to `kofamscan_results/`. For each `*_kofamscan.tsv`, produces two files:

| Output | Description |
|--------|-------------|
| `*.clean.tsv` | All predictions. Columns: `hc`, `gene_name`, `KO`, `threshold`, `score`, `e_value`, `ko_definition` |
| `*.clean.hc.tsv` | High-confidence only (score above KO threshold). Same columns minus `hc`. |

The `hc` column is `y`/`n` — `y` means the score exceeds the KO's predefined threshold (marked with `*` in KofamScan raw output).

Single-file mode:

```bash
./kofamscan_to_tsv.sh input.tsv output.tsv       # writes output.tsv and output.hc.tsv
./kofamscan_to_tsv.sh input.tsv                   # prints to stdout (clean TSV only)
```

### 4. Map KO to KEGG pathways

```bash
./kofamscan2pathway.sh
```

Downloads KO→pathway and pathway name tables from the [KEGG REST API](https://rest.kegg.jp), then maps both the full clean predictions and the high-confidence subset to KEGG pathways. Reference files are cached after the first download.

For each sample, produces:

| Output | Description |
|--------|-------------|
| `*.mapper.txt` | `gene_id<TAB>KO` from hc (no header) — ready for [KEGG Mapper Reconstruct](https://www.kegg.jp/kegg/mapper/reconstruct.html) web upload |
| `*.gene2pathway.tsv` | All predictions mapped to pathways: `gene_name<TAB>KO<TAB>pathway_id<TAB>pathway_name` |
| `*.gene2pathway.hc.tsv` | High-confidence only mapped to pathways (same columns) |

Example row:

```
HORVU.MOREX.r3.6HG0546330.1	K01951	map00230	Purine metabolism
```

Custom results directory:

```bash
./kofamscan2pathway.sh /path/to/results/
```

### KEGG Mapper web upload

The `.mapper.txt` files can be pasted directly into [KEGG Mapper Reconstruct](https://www.kegg.jp/kegg/mapper/reconstruct.html) to visualize which pathways are represented in your data.

## Output directory structure

```
kofamscan_results/
├── {sample}_kofamscan.tsv              # raw KofamScan output
├── {sample}_kofamscan.clean.tsv        # all predictions, clean headers
├── {sample}_kofamscan.clean.hc.tsv     # high-confidence only
├── {sample}_kofamscan.mapper.txt       # gene→KO (hc, for KEGG Mapper upload)
├── {sample}_kofamscan.gene2pathway.tsv    # all predictions → pathways
└── {sample}_kofamscan.gene2pathway.hc.tsv # high-confidence → pathways
```

## Reference files

Downloaded by `kofamscan2pathway.sh` and cached in the project root:

| File | Source |
|------|--------|
| `ko_pathway_link.tsv` | `https://rest.kegg.jp/link/pathway/ko` |
| `pathway_names.tsv` | `https://rest.kegg.jp/list/pathway` |

## Recreate conda environment

```bash
conda env create -f environment.yml
```
