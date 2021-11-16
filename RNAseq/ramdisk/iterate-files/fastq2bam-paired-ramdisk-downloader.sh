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

for sample in $(cat ${dwnlds} | tr -d '[]"' | tr , '\n')
do
   name= basename -- "${sample}" 
   p1=${sample}
   p2=${sample%%_1.fastq.gz}_2.fastq.gz
   wget $p1
   wget $p2
   bash fastq2bam-paired-ramdisk-iterative.sh -f ${name} -g Zm-B73-REFERENCE-NAM-5.0 -h Zmays_493_APGv4_Phytozome
done
