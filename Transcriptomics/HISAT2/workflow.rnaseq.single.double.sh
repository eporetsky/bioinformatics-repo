# $1 = experiment accession
# $2 = sample name
# $3 = fastq file 1
# $4 = genome index

axel --quiet -n 8 -o ramdisk $3

fastp --thread 15 --in1 ramdisk/$2.fastq.gz --out1 ramdisk/$2.fastq --html $1/reports/fastp_$2.html

hisat2 -p 15 --max-intronlen 6000 -x $4 -U ramdisk/$2.fastq --summary-file $1/reports/hisat2_$2.txt | \
sambamba view -S -f bam -o /dev/stdout /dev/stdin | \
sambamba sort -F "not unmapped" --tmpdir="ramdisk/tmpmba" -t 15 -o ramdisk/$2.bam /dev/stdin

featureCounts -t exon,CDS -T 15 -a $4.gtf -o $1/counts/$2.counts ramdisk/$2.bam

samtools view -@ 15 -T $4.fa -C -o $1/crams/$2.cram ramdisk/$2.bam

rm ramdisk/$2*
