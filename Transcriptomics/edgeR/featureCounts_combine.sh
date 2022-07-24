# Create a file containing the first column with the gene IDs. Skip 1st row which is a comment
# I used the first .counts file in the folder to get the GeneID column as the initial file
awk -F '\t' '(NR>1)' $(ls -tr *.counts | head -1) | awk -F '\t' '{print $1}' - > combined_counts.tsv

for name in *.counts
do
    # Combine the latest test.csv with the count column (last column $NF) of the next file 
    awk -F '\t' '(NR>1)' ${name} | awk -F '\t' '{print $NF}' - | paste -d'\t' combined_counts.tsv - > temp.tsv
    # Re-write the test.tsv file with the count data from the last file added
    cp temp.tsv combined_counts.tsv
done

rm temp.tsv
