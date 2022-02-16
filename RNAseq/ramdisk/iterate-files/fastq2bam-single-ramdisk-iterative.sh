while getopts f:g:b:h:t: flag

# This function downloads and processes one single-end fastq file at a time
# This is more useful for larger projects that will not fit on RAM without compression 
# TODO: Re-write it as a bash function in bash fastq2bam-paired-ramdisk-downloader.sh
# TODO: Integrate the md5checksum with the current function

do
    case "${flag}" in
    	f) filenm=${OPTARG};;
        g) gnzref=${OPTARG};;
        b) gnzbam=${OPTARG};;
        h) hstref=${OPTARG};;
        t) gtfref=${OPTARG};;
    esac
done

sample=${filenm%%.fastq.gz}
out=mapped/${sample}.bam.genozip

echo =========================================
echo Sample $sample
echo Genozip ref ${gnzref}.ref.genozip 
echo hisat2 ref ${hstref}.fa
echo =========================================
	
fastp --in1 ${sample}.fastq.gz --out1 "single.fastq.gz" --html reports/${sample}.html 
hisat2 -p 32 --max-intronlen 6000 -x ${hstref} -U single.fastq.gz |				          
sambamba view -S -f bam -o /dev/stdout /dev/stdin |
sambamba sort  --tmpdir="tmpmba" -t 32 -o mapped/${sample}.bam /dev/stdin

rm single.fastq.gz
rm ${sample}.fastq.gz

featureCounts -t exon,CDS -T 32 -a ${gtfref} -o counts/${sample}.counts mapped/${sample}.bam
genozip -e ${gnzbam} --replace -i bam -o $out mapped/${sample}.bam
rm mapped/${sample}.bam.bai


