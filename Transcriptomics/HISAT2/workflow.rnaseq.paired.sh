# $1 = experiment accession
# $2 = sample name
# $3 = fastq file 1
# $4 = md5 file 1
# $5 = fastq file 2
# $6 = md5 file 2
# $7 = genome index

echo fasterq-dump $2 --threads 36 -O ramdisk/
fasterq-dump $2 --threads 36 -O ramdisk/
mv ramdisk/$2_1.fastq ramdisk/tmp.$2_1.fastq
mv ramdisk/$2_2.fastq ramdisk/tmp.$2_2.fastq

fastp --thread 16 --in1 ramdisk/tmp.$2_1.fastq --out1 ramdisk/$2_1.fastq --in2 ramdisk/tmp.$2_2.fastq --out2 ramdisk/$2_2.fastq --html $1/reports/fastp_$2.html
rm ramdisk/tmp.$2_1.fastq ramdisk/tmp.$2_2.fastq

hisat2 -p 70 --min-intronlen 20 --max-intronlen 6000 -x $7 -1 ramdisk/$2_1.fastq -2 ramdisk/$2_2.fastq --summary-file $1/reports/hisat2_$2.txt | \
sambamba view -S -f bam -o /dev/stdout /dev/stdin | \
sambamba sort -F "not unmapped" --tmpdir="ramdisk/tmpmba" -t 72 -o ramdisk/$2.bam /dev/stdin

featureCounts -p --countReadPairs -t exon,CDS -T 35 -a $7.gtf -o $1/counts/$2.counts ramdisk/$2.bam

samtools view -@ 72 -T $7.fa -C -o $1/crams/$2.cram ramdisk/$2.bam

rm ramdisk/$2*
