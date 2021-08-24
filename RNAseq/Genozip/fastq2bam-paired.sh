#!/bin/bash

# Copied from: https://genozip.readthedocs.io/fastq-to-bam-pipeline.html
# TODO: Make notes and modifications to work with specified files

ref=GRCh38_full_analysis_set_plus_decoy_hla.fa
study=mystudy
fastq=myfastqdir
mapped=mymappeddir

files=($fastq/*.genozip)

processed=1

while (( processed == 1 )); do
    processsed=0
    for file in ${files[@]}
    do
        sample=`echo $file |grep -o -E sample'[[:digit:]]{2}'`  # convert the file name to a sample name
        out=$mapped/${sample}.bam.genozip

        if [ -f $out ]; then continue; fi # already processed
        if [ -f ${out}.doing_now ]; then continue; fi # another instance of this script is working on it

        processed=1
        touch ${out}.doing_now

        echo =========================================
        echo Sample $sample
        echo =========================================

        ( genocat --interleave $file -e ${ref%.fa}.ref.genozip                                                  || >&2 echo "genocat exit=$?" )|\
        ( fastp --stdin --interleaved_in --stdout --html ${fastq}/${sample}.html --json ${fastq}/${sample}.json || >&2 echo "fastp exit=$?"   )|\
        ( bwa mem $ref - -p -t 54 -T 0 -R "@RG\tID:$sample\tSM:$study\tPL:Illumina"                             || >&2 echo "bwa exit=$?"     )|\
        ( samtools view -h -OSAM                                                                                || >&2 echo "samtools exit=$?")|\
        ( bamsort fixmates=1 adddupmarksupport=1 inputformat=sam outputformat=sam inputthreads=5 outputthreads=5 sortthreads=30 level=1  || >&2 echo "bamsort exit=$?" )|\
        ( bamstreamingmarkduplicates inputformat=sam inputthreads=3 outputthreads=3 level=1                     || >&2 echo "bamstreamingmarkduplicates exit=$?" )|\
        ( genozip -e ${ref%.fa}.ref.genozip -i bam -o $out -t                                                   || >&2 echo "genozip exit=$?" )

        rm ${out}.doing_now
    done
done