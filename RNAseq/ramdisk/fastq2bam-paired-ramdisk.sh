#!/bin/bash
# Modified from: https://genozip.readthedocs.io/fastq-to-bam-pipeline.html
# Bash scripts that works with paired-end data. Produces bam.genozip for each sample.

while getopts g:h: flag
do
    case "${flag}" in
        g) gnzref=${OPTARG};;
        h) hstref=${OPTARG};;
    esac
done

mkdir mapped/
mkdir reports/
mkdir counts/

for file in *_1.fastq.gz
do
    sample=${file%%_1.fastq.gz}
    out=mapped/${sample}.bam.genozip

    echo =========================================
    echo Sample $sample
    echo File $file
    echo Genozip ref ${gnzref}.ref.genozip 
    echo hisat2 ref ${hstref}.fa
    echo =========================================
	
    #genocat --interleave $file -e ${gnzref}.ref.genozip |
    fastp --in1 ${sample}_1.fastq.gz --out1 "pair1.fastq.gz" --in2 ${sample}_2.fastq.gz --out2 "pair2.fastq.gz" --html reports/${sample}.html 
    hisat2 -p 32 --max-intronlen 6000 -x ${hstref}.fa -1 pair1.fastq.gz -2 pair2.fastq.gz |				          
    sambamba view -S -f bam -o /dev/stdout /dev/stdin |
    sambamba sort  --tmpdir="tmpmba" -t 32 -o mapped/${sample}.bam /dev/stdin
    featureCounts -t exon,CDS -T 32 -a Zmays_493_RefGen_V4_Phytozome.gene.gtf -o counts/${sample}.counts mapped/${sample}.bam
    genozip -e ${hstref}.ref.genozip --replace -i bam -o $out mapped/${sample}.bam
    rm mapped/${sample}.bam.bai
    genozip --replace --reference ${gnzref}.ref.genozip --pair ${sample}_1.fastq.gz ${sample}_2.fastq.gz
done

rm pair1.fastq.gz
rm pair2.fastq.gz
