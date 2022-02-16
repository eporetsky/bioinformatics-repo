while getopts d:g:b:h:t: flag
do
    case "${flag}" in
        d) dwnlds=${OPTARG};;
        g) gnzref=${OPTARG};;
        b) gnzbam=${OPTARG};;
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
   axel -a -n 32 ${sample}
   bash fastq2bam-single-ramdisk-iterative.sh -f ${name} -g ${gnzref} -b ${gnzbam} -h ${hstref} -t ${gtfref}
done
