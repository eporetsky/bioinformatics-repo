while getopts r:g: flag
do
    case "${flag}" in
        r) ref=${OPTARG};;
        g) gtf=${OPTARG};;
    esac
done

for file in bams/*.bam.genozip
do
    bam=${file%%.genozip}
    
    sample=$(basename ${file}) 
    sample=${sample%%.bam.genozip}
    out=counts/${sample}.counts

    echo =========================================
    echo Sample $sample
    echo File $file
    echo Genozip ref ${ref} 
    echo GTF file ${gtf}
    echo =========================================
    
    genounzip --replace --reference ${ref} ${file}
    featureCounts  -p -t exon,CDS -T 32 -a ${gtf} -o ${out} ${bam}
    rm ${bam}
done
