# $1 = experiment accession
# $2 = sample name
# $3 = genome index

# Function to calculate the MD5 checksum of the downloaded file

fastp --thread 16 --in1 $1/fastq/$2_1.fq.gz --out1 ramdisk/$2_1.fastq --in2 $1/fastq/$2_2.fq.gz --out2 ramdisk/$2_2.fastq --html $1/reports/fastp_$2.html

hisat2 -p 70 --min-intronlen 20 --max-intronlen 6000 -x $3 -1 ramdisk/$2_1.fastq -2 ramdisk/$2_2.fastq --summary-file $1/reports/hisat2_$2.txt | \
sambamba view -S -f bam -o /dev/stdout /dev/stdin | \
sambamba sort -F "not unmapped" --tmpdir="ramdisk/tmpmba" -t 35 -o ramdisk/$2.bam /dev/stdin

featureCounts -p --countReadPairs -t exon,CDS -T 35 -a $3.gtf -o $1/counts/$2.counts ramdisk/$2.bam

samtools view -@ 36 -T $3.fa -C -o $1/crams/$2.cram ramdisk/$2.bam

rm ramdisk/$2*