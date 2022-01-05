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

## Analyze one sample at a time using the bash script in the iterate-files
```
# For example:
-d: txt file with fastq ftp links, one line per file. Should add ftp:// to beginning of line of missing
-g: Genozip reference file. I use the same one for all fastq files of the same species. 
-h: Name of hisat2 index and alignment-specific genozip file, use same genome for both.
-t: GTF file that will be used with feature counts.
bash fastq2bam-single-ramdisk-downloader.sh -d ftp_file_list.txt -g Zm-B73-REFERENCE-NAM-5.0 -h Zmays_493_APGv4_Phytozome -t Zmays_493_RefGen_V4_Phytozome.gene.gtf
```

## Alternatively: Download all the files to the ramdisk
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
