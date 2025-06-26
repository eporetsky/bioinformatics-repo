# Docking Utilities

This subdirectory contains scripts and utilities for protein-ligand docking workflows.

## Installation

To set up the required environment, use the provided `environment.yml` file. This environment is based on the [LiganDock](https://github.com/eporetsky/ligandock) repository and includes all necessary dependencies for the scripts in this directory.

```bash
conda env create -f environment.yml
conda activate ligandock
```
---

## Script: `combine.ligands.py`

### Description

`combine.ligands.py` is a utility script for merging multiple ligand SDF files and/or metal PDB files into a single SDF file. This is useful for preparing input files for docking workflows where you need to combine organic ligands and metal ions into one structure.

- **Ligand SDF files** are used as-is.
- **Metal PDB files** are converted to SDF format using RDKit before merging.
- The script requires at least two input molecules (from SDF and/or PDB files).

### Usage

```bash
python combine.ligands.py --input_sdf ligand1.sdf ligand2.sdf --input_pdb metal1.pdb metal2.pdb --output_sdf merged_output.sdf
```

#### Arguments
- `--input_sdf`: One or more input ligand SDF files (optional, can be empty if only PDBs are provided)
- `--input_pdb`: One or more input metal PDB files (optional, can be empty if only SDFs are provided)
- `--output_sdf`: Output file name for the merged SDF (required)

#### Example
Combine two ligands and one metal ion:
```bash
python combine.ligands.py \
    --input_sdf ligand1.sdf ligand2.sdf \
    --input_pdb zn.pdb --output_sdf merged.sdf
```

---