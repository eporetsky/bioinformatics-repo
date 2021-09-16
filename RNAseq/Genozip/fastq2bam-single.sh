#!/bin/bash
# Modified from: https://genozip.readthedocs.io/fastq-to-bam-pipeline.html
# Bash scripts that works with single-end data. Produces bam.genozip for each sample.

while getopts g:b: flag
do
    case "${flag}" in
        g) gnzref=${OPTARG};;
        b) bwaref=${OPTARG};;
    esac
done

#gnzref=Zm-B73-REFERENCE-NAM-5.0 
#bwaref=Zmays_493_APGv4_Phytozome

mkdir mapped/
mkdir reports/

for file in *.genozip
do
    sample=`echo $file | awk -F'[_.]' '{print $1}'` # if the files are in a different folder: | awk -F'[_/]' '{print $2}'`
    out=mapped/${sample}.bam.genozip

    echo =========================================
    echo Sample $sample
    echo File $file
    echo Genozip ref ${gnzref}.ref.genozip 
    echo Bwa-mem ref ${bwaref}.fa
    echo =========================================

    genocat $file -e ${gnzref}.ref.genozip |
    fastp --stdin --stdout --stdout --html reports/${sample}.html  |
    bwa mem ${bwaref}.fa - -t 32 |
    sambamba view -S -f bam -o /dev/stdout /dev/stdin |
    sambamba sort  --tmpdir="tmpmba" -t 32 -o /dev/stdout /dev/stdin |
    genozip -e ${bwaref}.ref.genozip -i bam -o $out -t                          

done

