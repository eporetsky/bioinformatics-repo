while getopts f:g:h: flag

# This function downloads and processes one pair of fastq files at a time
# This is more useful for larger projects that will not fit on RAM without compression 
# TODO: Re-write it as a bash function in bash fastq2bam-paired-ramdisk-downloader.sh
# TODO: Integrate the md5checksum with the current function

do
    case "${flag}" in
    	f) filenm=${OPTARG};;
        g) gnzref=${OPTARG};;
        h) hstref=${OPTARG};;
    esac
done

sample=${filenm%%_1.fastq.gz}
out=mapped/${sample}.bam.genozip

echo =========================================
echo Sample $sample
echo Genozip ref ${gnzref}.ref.genozip 
echo hisat2 ref ${hstref}.fa
echo =========================================
	
fastp --in1 ${sample}_1.fastq.gz --out1 "pair1.fastq.gz" --in2 ${sample}_2.fastq.gz --out2 "pair2.fastq.gz" --html reports/${sample}.html 
hisat2 -p 32 --max-intronlen 6000 -x ${hstref}.fa -1 pair1.fastq.gz -2 pair2.fastq.gz |				          
sambamba view -S -f bam -o /dev/stdout /dev/stdin |
sambamba sort  --tmpdir="tmpmba" -t 32 -o mapped/${sample}.bam /dev/stdin

rm pair1.fastq.gz
rm pair2.fastq.gz

featureCounts -p -t exon,CDS -T 32 -a Zmays_493_RefGen_V4_Phytozome.gene.gtf -o counts/${sample}.counts mapped/${sample}.bam
genozip -e ${hstref}.ref.genozip --replace -i bam -o $out mapped/${sample}.bam
rm mapped/${sample}.bam.bai
genozip --replace --reference ${gnzref}.ref.genozip --pair ${sample}_1.fastq.gz ${sample}_2.fastq.gz


