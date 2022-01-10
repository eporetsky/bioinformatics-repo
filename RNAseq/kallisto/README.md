# Analysis of 3' RNA-seq results using Kallisto
The following code is adapted from the tutorial:\
https://pachterlab.github.io/kallisto/starting\
\
While this is not stricktly required, I wanted to ensure the presense of a minimal 500bp extension after the stop-codon.\
The script `fast3p_extender.ipynb` can generated the new cDNA files that will be used to generate the Kallisto index.\
The script uses the genome FASTA and GFF file to generate a single sequence for each gene made of the last 500bp of the\
CDS sequence + the annotated 3' UTR (if exists) + extension to 500bp following the stop codon (in required).

## Required files 
* Kallisto (conda install -c bioconda kallisto)
* Genome FASTA
* Genome GFF3 
* BioPython

## To generate a Kallisto index from a reference cDNA fasta file
```
for d in *.fa ; do (kallisto index -i ${d%_fast3p.fa}.idx $d) ; done
```

## FASTQ filtering using fastp
```
mkdir filtered
mkdir reports
for d in *.fastq.gz ; do (fastp --in1 $d --out1 filtered/${d%.fastq.gz}.fp.fastq.gz --html reports/${d%.fastq.gz}.html --thread 16 --compression 7) ; done
```

## To download and compress files using genozip:
```
wget -i ftp_file_list.txt
```

## Run the Kallisto analysis
Python code for quantification using Kallisto and data parsers see: https://github.com/eporetsky/SQNce/blob/main/helpers/TPM_tables.ipynb

## Notes
* Work with pseudobam files: https://www.biostars.org/p/263589/ (alternatively treat every gene as chromosome)
* Check existing parsers of Kallisto output: https://rdrr.io/bioc/SummarizedExperiment/man/readKallisto.html (I wrote my own parser for now)
* Include bootstrap? https://nbisweden.github.io/workshop-RNAseq/1811/labs/kallisto.html (Bootstraps is mandatory for isoforms differential expression analysis with Sleuth)
* TODO: Include a link to dataset quantification
