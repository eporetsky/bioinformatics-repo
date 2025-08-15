# RNA-seq Analysis Workflow

This project provides a workflow for processing RNA-seq data using publicly available datasets from the European Nucleotide Archive (ENA). The analysis is designed to download and process experiments based on directory names that match ENA accessions, such as PRJNA1070606.

## Installation

To set up the environment, use the provided `environment.yml` file to create a conda environment:

```bash
conda env create --file environment.yml
conda activate rnaseq
```

## Optional: Creating a RAM Disk

For users with more than 100GB of RAM, creating a RAM disk can significantly speed up the analysis by loading files into RAM instead of the hard disk, while sparing the hard drive from excessive writing. This step is optional and should only be performed if you have sufficient RAM.

```bash
mkdir ramdisk
sudo mount -t tmpfs -o size=100000m tmpfs ramdisk/
```

## Running the Analysis

1. **Prepare the Directory Structure**: Ensure that each experiment is in its own directory, named after the ENA accession number.
2. **Download and Process Data**: The scripts will automatically download and process the data based on the directory names.

## Data Requirements

- Directory names should match ENA accession numbers.
- Ensure that the necessary HISAT2 genome index and GTF files are available in the working directory.
- The HISAT2 and GTF files should have the name for the workflow to work.
- The workflow handles experiments with mixed single- and paired-end FASTQ files individually
- In samples with both single- and paired-end FASTQ files, the single-end files are ignored

## Scripts Overview

- `workflow.rnaseq.paired.sh`: Processes paired-end RNA-seq data.
- `workflow.rnaseq.single.sh`: Processes single-end RNA-seq data.
- `workflow.rnaseq.py`: Manages the overall workflow, downloading data and executing the appropriate scripts.
- `combine.counts.sh`: Combines count data from multiple experiments into a single file.

## Example Usage

To run the analysis for a specific experiment:

```bash
# Process alldirectories with proper accession IDs
python workflow.rnaseq.py

# Process all the individual count files into one file
bash combine.counts.sh 
```

## Troubleshooting

- Ensure that the conda environment is activated before running the scripts.
- Check that the RAM disk is properly mounted if using that option.
- Verify that directory names match ENA accession numbers and that all required files are present.
- Supposedely, fasterq-dump runs checks on each downloaded chunk
- Users should manually check all samples were downloaded correctly 