#!/bin/bash

# Modified from: https://genozip.readthedocs.io/fastq-to-bam-pipeline.html
# Bash scripts that works with single-end data. Produces bam.genozip for each sample.

ref=reference_genome_name   # Assume suffix is .ref.genozip (do not include in name)
files=(genozip/*.genozip)   # Each fastq file is compressed as a genozip file

processed=1

while (( processed == 1 )); do
    processsed=0
    for file in genozip/*.genozip
    do
        # Convert the file name to a sample name based on the below assumptions
        # Assumes no `.` in folder names and assumed existance of `/` in $file
        sample=`echo $file | awk -F'[_.]' '{print $1}' | awk -F'[_/]' '{print $2}'`
        out=mapped/${sample}.bam.genozip

        if [ -f $out ]; then continue; fi # already processed
        if [ -f ${out}.doing_now ]; then continue; fi # another instance of this script is working on it

        processed=1
        touch ${out}.doing_now

        echo =========================================
        echo Sample $sample
        echo =========================================

        ( genocat $file -e ${ref}.ref.genozip                                                  	         	    || >&2 echo "genocat exit=$?" )|\
        ( fastp --stdin --stdout --html genozip/${sample}.html || >&2 echo "fastp exit=$?"   )|\
        ( bwa mem reference_genome_name.bw.ix - -t 54 -T 0                             			      		    || >&2 echo "bwa exit=$?"     )|\
        ( samtools view -h -OSAM                                                                                || >&2 echo "samtools exit=$?")|\
        ( bamsort fixmates=1 adddupmarksupport=1 inputformat=sam outputformat=sam inputthreads=5 outputthreads=5 sortthreads=30 level=1  || >&2 echo "bamsort exit=$?" )|\
        ( bamstreamingmarkduplicates inputformat=sam inputthreads=3 outputthreads=3 level=1                     || >&2 echo "bamstreamingmarkduplicates exit=$?" )|\
        ( genozip -e ${ref%.fa}.ref.genozip -i bam -o $out -t                                                   || >&2 echo "genozip exit=$?" )

        rm ${out}.doing_now
    done
    # Script currently gets stuck in an infinite loop when done. TODO: Fix break points.
done
