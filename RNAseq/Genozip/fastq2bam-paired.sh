#!/bin/bash

ref="genome_name"

mkdir mapped/
mkdir reports/

for file in *.genozip
do
    sample=`echo $file | awk -F'[_.]' '{print $1}'` # if the files are in a different folder: | awk -F'[_/]' '{print $2}'`
    out=mapped/${sample}.bam.genozip

    echo =========================================
    echo Sample $sample
    echo File $file
    echo Ref ${ref}.ref.genozip 
    echo =========================================

    genocat --interleave $file -e ${ref}.ref.genozip                       |
    fastp --stdin --stdout --html --interleaved_in --stdout --html reports/${sample}.html  |
    bwa mem $ref - -p -t 32 |					          
    sambamba view -S -f bam -o /dev/stdout /dev/stdin |
    sambamba sort  --tmpdir="tmpmba" -t 32 -o /dev/stdout /dev/stdin |
    genozip -e ${ref}.ref.genozip -i bam -o $out -t                          

    done

