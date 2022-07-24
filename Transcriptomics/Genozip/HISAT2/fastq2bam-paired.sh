#!/bin/bash
# Modified from: https://genozip.readthedocs.io/fastq-to-bam-pipeline.html
# Bash scripts that works with paired-end data. Produces bam.genozip for each sample.



while getopts g:b: flag
do
    case "${flag}" in
        g) gnzref=${OPTARG};;
        h) hstref=${OPTARG};;
    esac
done

mkdir mapped/
mkdir reports/

for file in *.genozip
do
    sample=${file%%_1+2.fastq.genozip}
    out=mapped/${sample}.bam.genozip

    echo =========================================
    echo Sample $sample
    echo File $file
    echo Genozip ref ${gnzref}.ref.genozip 
    echo hisat2 ref ${hstref}.fa
    echo =========================================
	
    genocat --interleave $file -e ${gnzref}.ref.genozip                      |
    fastp --stdin --out1 "pair1.fastq.gz" --out2 "pair2.fastq.gz" --interleaved_in --html reports/${sample}.html 
    hisat2 -p 32 --max-intronlen 6000 -x ${hstref}.fa -1 pair1.fastq.gz -2 pair2.fastq.gz |				          
    sambamba view -S -f bam -o /dev/stdout /dev/stdin |
    sambamba sort  --tmpdir="tmpmba" -t 32 -o /dev/stdout /dev/stdin |
    genozip -e ${hstref}.ref.genozip -i bam -o $out                        

done

rm pair1.fastq.gz
rm pair2.fastq.gz
