#!/bin/bash

rm -r test/ 2> /dev/null
rm -r pMT/dist/ 2> /dev/null


awk -F '\t' 'NR==FNR {id[$1]; next} $1 in id' AtIDs.txt Compara.94.protein_default.homologies.tsv > compara_temp.tsv 
awk -F '\t' 'NR==FNR {id[$1]; next} $8 in id' PlantIDs.txt compara_temp.tsv >compara_final.tsv
awk '{print $2"\t"$8"\t"$6"\t"$4}' compara_final.tsv | sort  -k1 -k2 -rk4 | sort -u -k1,2 > top_orthologs.tsv

mkdir test 
awk '{print "test/"$1}' top_orthologs.tsv | uniq | xargs -I % touch %.fasta

cd test 
find *.fasta | 
sed -e 's/\.[^.]*$//' | 
xargs -P 20 -I % sh -c 'curl "http://rest.ensemblgenomes.org/sequence/id/%?type=protein" -H "Content-type:text/x-fasta" | 
sed "1 s/.*/&_arabidopsis_thaliana /"  >> %.fasta'
cd ..

# Why am I running this instead of on the top_orthologs.tsv
#awk '{print $2"\t"$8"\t"$7"\t"$4"\t"$8}' compara_final.tsv |
awk '{print $2"\t"$8"\t"$7"\t"$4"\t"$8}' top_orthologs.tsv | 
sort  -k1 -k2 -rk4 | 
sort -u -k1,2 | 
awk '{print $1"\t"$3"\t"$5}' | 
xargs -P 50 -I % sh -c 'echo % | 
(read At Orth OrthSpecies ; 
curl "http://rest.ensemblgenomes.org/sequence/id/$Orth?type=protein" -H 'Content-type:text/x-fasta' | 
tr "_" "-" |
sed "1 s/.*/&_\\"$OrthSpecies"/" >> test/"$At".fasta)'

cd test
find *.fasta | 
sed -e 's/\.[^.]*$//' | 
xargs -P 20 -I % sh -c 'do ./../muscle3.8.31_i86darwin64.fasta -in $.fast -out %.aln'
cd ..

cd test
find *.fasta | 
sed -e 's/\.[^.]*$//' | 
xargs -P 20 -I % sh -c './../muscle3.8.31_i86darwin64.fasta -in %.fasta -out %.aln'

find *.aln | 
sed -e 's/\.[^.]*$//' | 
xargs -P 20 -I % sh -c 'treebest nj %.aln> %.nhx'

find *.aln | 
sed -e 's/\.[^.]*$//' | 
xargs -P 20 -I % sh -c 'treebest distmat jtt %.aln | awk "{if(NR>1)print}" > %.tempdist'
cd ..

mkdir pMT/dist/
cp test/*.tempdist pMT/dist/

Rscript pMT/tempdist2dist.R

rm pMT/dist/*.tempdist

#awk '{print $8}' compara_temp.tsv | uniq >> pMT/alltaxids.txt
Rscript pMT/mt_pvals.R pMT/alltaxids.txt pMT/dist/ 10 0.05 1000 results.csv


#awk '{ if (($5 < 0.001) && ($3 > 20)) { print } }' results.csv > results_filtered.csv
#awk '{ if (($5 < 0.05)) { print } }' results.csv > results_filtered.csv
#cut -f1 results_filtered.csv >> filter_degrees.txt
#cut -f2 results_filtered.csv >> filter_degrees.txt
#sort filter_degrees.txt | uniq -c | awk '(($1<=10) && ($1>1)){print $2}' > results_keep.txt
#'NR==FNR{F1[$0];next}$1 in F1{print}' results_keep.txt results_filtered.csv >> level10.csv
#'NR==FNR{F1[$0];next}$2 in F1{print}' results_keep.txt results_filtered.csv >> level10.csv

#for i in *.fastq.gz
#do
#  input=$i
#  output="bams/"${i%%.fastq.gz}".bam"
#  hisat2 -x maizeV4_index -p 4 --min-intronlen 60 --max-intronlen 6000 -U ${input} |samtools view -bS > ${output}
#done
