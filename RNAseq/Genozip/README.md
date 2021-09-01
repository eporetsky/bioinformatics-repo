# RNA-seq Analysis using Genozip compressed files
Use the included files to convert fastq to bam files\
Adapted from: https://genozip.readthedocs.io/fastq-to-bam-pipeline.html\

## Required files 
TODO: Complete section
* Genozip
* bwa
* quant3p
* snakemake
* few more missing

## To download and compress files using genozip:
```
wget -i ftp_file_list.txt
# paired-end
for d in *1.fastq.gz ; do ( genozip --replace --reference genome.ref.genozip --pair $d ${d%_1.fastq.gz}_2.fastq.gz); done
# single-end (not tested yet, but should work)
for d in *.fastq.gz ; do ( genozip --replace --reference genome.ref.genozip $d; done
```

## To quantify the bam.genozip files using quant3p:
NOTE: Used for 3' RNA-seq data only - https://github.com/ctlab/quant3p\
TODO: Confirm that the code below works\
\
To run rTASSEL GLM with GNU parallel:\
(The goal of using ramdisk is to reduce writing to SSD drive+speed-up analysis)\
Make ramdisk in the folder containing the bam files
```
mkdir ramdisk
```
Mount ramdisk of specified size
```
sudo mount -t tmpfs -o size=30000m tmpfs ramdisk/
```
Unzip and create index using genozip
```
for d in find *.genozip ; do (genounzip $d --index --reference reference.ref.genozip --output ramdisk/${d%.genozip}); done
```
Run the quant3p analysis pipeline
```
cd ramdisk
quant3p -p 32 -n project_name -g reference.gtf *.bam
```

## Notes
* bwa mem doesn't recognized paired-end reads if genozip with --optimize-DESC
