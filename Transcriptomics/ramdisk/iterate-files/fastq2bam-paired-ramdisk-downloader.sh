while getopts d:g:h: flag
do
    case "${flag}" in
        d) dwnlds=${OPTARG};;
        g) gnzref=${OPTARG};;
        h) hstref=${OPTARG};;
    esac
done

mkdir mapped/
mkdir reports/
mkdir counts/

# Iterate over the downloads.txt file that contains the file names
for sample in $(cat ${dwnlds} | tr -d '[]"' | tr , '\n')
do
   name=$(basename ${sample}) 
   p1=${sample}
   p2=${sample%%_1.fastq.gz}_2.fastq.gz
   wget $p1
   wget $p2
   bash fastq2bam-paired-ramdisk-iterative.sh -f ${name} -g ${gnzref} -h ${hstref}
done
