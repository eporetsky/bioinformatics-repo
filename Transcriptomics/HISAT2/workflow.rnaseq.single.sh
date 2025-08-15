# $1 = experiment accession
# $2 = sample name
# $3 = fastq file 1
# $4 = md5 file 1
# $5 = genome index

fasterq-dump $2 --threads 36 -O ramdisk/
mv ramdisk/$2.fastq ramdisk/tmp.$2.fastq

fastp --thread 16 --in1 ramdisk/tmp.$2.fastq --out1 ramdisk/$2.fastq --html $1/reports/fastp_$2.html
rm ramdisk/tmp.$2.fastq

hisat2 -p 70 --max-intronlen 6000 -x $5 -U ramdisk/$2.fastq --summary-file $1/reports/hisat2_$2.txt | \
sambamba view -S -f bam -o /dev/stdout /dev/stdin | \
sambamba sort -F "not unmapped" --tmpdir="ramdisk/tmpmba" -t 72 -o ramdisk/$2.bam /dev/stdin

featureCounts -t exon,CDS -T 35 -a $5.gtf -o $1/counts/$2.counts ramdisk/$2.bam

samtools view -@ 72 -T $5.fa -C -o $1/crams/$2.cram ramdisk/$2.bam

rm ramdisk/$2*
