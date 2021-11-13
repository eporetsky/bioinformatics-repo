# Started using ramdisk on linux to download and process RNA-seq raw data
Use the included files to convert fastq to bam files\
Adapted from: https://genozip.readthedocs.io/fastq-to-bam-pipeline.html\

## Required files 
TODO: Complete section
* Genozip
* hisat2
* fastp
* sambamba
* featureCounts

## Create the ramdisk
(The goal of using ramdisk is to reduce writing to SSD drive+speed-up analysis)\
Make ramdisk in the folder containing the bam files
```
cd mapped
mkdir ramdisk
```
Mount ramdisk of specified size
```
sudo mount -t tmpfs -o size=80000m tmpfs ramdisk/
```

## Download all the files to the ramdisk
```
wget -i ftp_file_list.txt
# paired-end
for d in *_1.fastq.gz ; do ( genozip --replace --reference genome.ref.genozip --pair $d ${d%_1.fastq.gz}_2.fastq.gz); done
# single-end
for d in *.fastq.gz ; do ( genozip --replace --reference genome.ref.genozip $d ) ; done
# # If your paired-end data comes with singleton fastq files, manually rename to *.siingleton.fastq.gz for later use
```

## Run the bashscript
```
# paired-end (not tested yet but should work)
bash fastq2bam-paired-ramdisk.sh -g genozip_ref_name -h hisat2_index_name
# single-end or singletons (in progress)
# bash fastq2bam-single.sh -g genozip_ref_name -h hisat2_index_name
```

## Notes
All the files are compressed with genozip after alignment and raw read counting\
Save files to hard-drive or NAS drive after analysis is complete to not lose it.
