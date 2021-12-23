while getopts d:g:h:t: flag
do
    case "${flag}" in
        d) dwnlds=${OPTARG};;
        g) gnzref=${OPTARG};;
        h) hstref=${OPTARG};;
        t) gtfref=${OPTARG};;
    esac
done

mkdir mapped/
mkdir reports/
mkdir counts/

# Iterate over the downloads.txt file that contains the file names
for sample in $(cat ${dwnlds} | tr -d '[]"' | tr , '\n')
do
   name=$(basename ${sample}) 
   wget ${sample}
   bash fastq2bam-paired-ramdisk-iterative.sh -f ${name} -g ${gnzref} -h ${hstref} -t ${gtfref}
done
