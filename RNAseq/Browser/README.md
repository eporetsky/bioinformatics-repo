# Prepare and visualize aligned RNA-seq data from bam.genozip files
Adapted from: https://genozip.readthedocs.io/\


## (The goal of using ramdisk is to reduce writing to SSD drive+speed-up analysis)\
Make ramdisk in the folder containing the bam files
```
cd mapped
mkdir ramdisk
# Mount ramdisk of specified size
sudo mount -t tmpfs -o size=80000m tmpfs ramdisk/
```

## To unzip the bam files and create an index file
```
for d in *.genozip ; do (genounzip $d --index --replace --reference genome.ref.genozip --output ${d%.genozip}); done
```

## To merge the files
The goal is to visualize the merged bam files from each condition
```
samtools merge merged.bam file1.bam file2.bam...
```
