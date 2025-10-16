#!/bin/bash
# ls clean > process.tsv

rm process.tsv

# Iterate over all .fa files in the input directory
for fasta_file in splits/*.fasta; do
    output_file="output/$(basename "$fasta_file").tsv"
    if [ ! -f "$output_file" ]; then
        echo $fasta_file >> process.tsv
    fi
done

mkdir output

parallel -j20 --colsep '\\t' -a process.tsv bash interproscan.sh -i {1} --output-dir output --formats TSV --iprlookup --goterms --disable-residue-annot --applications SUPERFAMILY,PANTHER,Gene3D,Pfam,FunFam,NCBIfam
# --pathways

#--applications CDD,Coils,FunFam,Gene3D,MobiDBLite,NCBIfam,PANTHER,Pfam,SUPERFAMILY
