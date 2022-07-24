#!/bin/bash

# bash script that aligns filtered fastq files (I use fastp, hence the name) using HISAT2 and generates BAM files using SAMBAMBA
# Use the -h flag to input the HISAT2 index name

while getopts h: flag
do
    case "${flag}" in
        h) hstref=${OPTARG};;
    esac
done

mkdir mapped/

for file in *.fastq.gz
do
    sample=${file%%.fastq.gz}
    out=mapped/${sample}.bam

    echo =========================================
    echo Sample $sample
    echo File $file
    echo hisat2 ref ${hstref}
    echo =========================================
	
    hisat2 -p 32 --max-intronlen 6000 -x ${hstref} -U $file |				          
    sambamba view -S -f bam -o /dev/stdout /dev/stdin |
    sambamba sort  --tmpdir="tmpmba" -t 32 -o $out /dev/stdin  
    sambamba index -t 32 $out           

done
