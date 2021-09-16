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

## To generate a genozip reference genome file
```
genozip --make-reference reference_genome.fa.gz
```

## To download and compress files using genozip:
```
wget -i ftp_file_list.txt
# paired-end
for d in *_1.fastq.gz ; do ( genozip --replace --reference genome.ref.genozip --pair $d ${d%_1.fastq.gz}_2.fastq.gz); done
# single-end
for d in *.fastq.gz ; do ( genozip --replace --reference genome.ref.genozip $d ) ; done
# # If your paired-end data comes with singleton fastq files, manually rename to *.siingleton.fastq.gz for later use
```

## To run the bwa mem aligner using genozip:
```
# paired-end (not tested yet but should work)
bash fastq2bam-paired.sh -g genozip_ref_name -b bwa_index_name
# single-end or singletons
bash fastq2bam-single.sh -g genozip_ref_name -b bwa_index_name

# genozip_ref_name refers to the original fastq genozip fastq compression reference
# bwa_index_name refers to the name of the current bwa index and genozip reference files
#   because bams can't be compressed with genozip.ref different to the alignment index 
```

## To quantify the bam.genozip files:
TODO: Confirm that the code below works\
\
Make sure you have a gtf file. To convert gff3 to gtf:
```
gffread reference.gff3 -T -o reference.gtf
```

(The goal of using ramdisk is to reduce writing to SSD drive+speed-up analysis)\
Make ramdisk in the folder containing the bam files
```
cd mapped
mkdir ramdisk
```
Mount ramdisk of specified size
```
sudo mount -t tmpfs -o size=30000m tmpfs ramdisk/
cd ramdisk
```
Unzip and create index using genozip
```
for d in find *.genozip ; do (genounzip $d --index --reference reference.ref.genozip --output ramdisk/${d%.genozip}); done

# If your paired-end data comes with singleton fastq files, merge the two files using samtools
cd ramdisk/
for d in *.singleton.bam ; do(samtools merge --write-index -o ${d%.singleton.bam}.merged.bam ${d%.singleton.bam}.bam $d); done
# Delete the original bam files so they are not included in the count. For loop deletes all files if no singleton.bam is present
for d in *.singleton.bam ; do(rm ${d%.singleton.bam}.bam ${d%.singleton.bam}.bam.bai $d $d.bai); done
```

#### To quantify the bam.genozip files using featureCounts:
```
featureCounts -T 32 -a ../../reference.gtf -o counts.txt *.bam
```


#### To quantify the bam.genozip files using quant3p:
NOTE: Used for 3' RNA-seq data only - https://github.com/ctlab/quant3p \
Run the quant3p analysis pipeline
```
cd ramdisk
quant3p -p 32 -n project_name -g reference.gtf *.bam
```

## Unmount ramdisk when done
Note: Don't forget to backup any files before unmounting or turning off computer

## Notes
* bwa mem doesn't recognized paired-end reads if genozip with --optimize-DESC
