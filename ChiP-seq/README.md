# ChiP-seq Analysis Workflow

I have been working on a ChIP-seq analysis tutorial on my website: [https://eporetsky.github.io/tutorials/chipseq](https://eporetsky.github.io/tutorials/chipseq). This page will contain information about specific scripts that I wrote for the tutorial. 

# diffbind2fasta.py

This a Python script designed to extract DNA sequences from a genome FASTA file based on genomic coordinates provided in a BED file. I assume it will work with any BED file but I have only tested it with the BED file generated with DiffBind, hence the name. This is especially useful when you just want to get the sequences from your own genome FASTA file and you don't want to use GRanges.

- Extract sequences from a genome FASTA file using genomic coordinates from a BED file.
- Automatically handles reverse complement sequences if the strand information is provided.
- Outputs the extracted sequences in FASTA format for downstream analysis.
- The script only requires BioPython

The script requires three input files:
1. A BED file (e.g., DiffBind output) containing genomic regions of interest.
2. A genome reference FASTA file (e.g., human genome or any organism of interest).
3. An output file name to save the extracted sequences in FASTA format.

Running the script:

```
python diffbind2fasta.py -b <path_to_bed_file> -f <path_to_genome_fasta> -o <path_to_output_fasta>
```

To reproduce the BED files I've created with DiffBind:

```
dba_report <- dba.report(dba_obj, method = DBA_EDGER, th = 0.05, bUsePval = FALSE, DataType = DBA_DATA_GRANGES)  
rtracklayer::export(dba_report, con = "output.bed", format = "BED")
```