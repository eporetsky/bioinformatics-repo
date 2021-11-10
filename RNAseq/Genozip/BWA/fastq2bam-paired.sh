#!/bin/bash
# Modified from: https://genozip.readthedocs.io/fastq-to-bam-pipeline.html
# Bash scripts that works with paired-end data. Produces bam.genozip for each sample.



while getopts g:b: flag
do
    case "${flag}" in
        g) gnzref=${OPTARG};;
        b) bwaref=${OPTARG};;
    esac
done

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

    genocat --interleave $file -e ${gnzref}.ref.genozip                      |
    fastp --stdin --stdout --html --interleaved_in --stdout --html reports/${sample}.html  |
    bwa mem ${bwaref}.fa - -p -t 32 |					          
    sambamba view -S -f bam -o /dev/stdout /dev/stdin |
    sambamba sort  --tmpdir="tmpmba" -t 32 -o /dev/stdout /dev/stdin |
    genozip -e ${bwaref}.ref.genozip -i bam -o $out -t                          

    done

