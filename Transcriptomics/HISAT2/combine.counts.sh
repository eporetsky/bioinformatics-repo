#!/bin/bash

# Iterate over each directory in the current directory
for experiment in */
do

    # Check if the directory name is "ramdisk/"
    if [[ "$experiment" == "ramdisk/" ]]; then
        echo "Skipping: $experiment"
        continue
    fi
    
    echo Starting: $experiment
    # Remove the trailing slash from the directory name
    experiment=${experiment%/}

    # Create a file containing the first column with the gene IDs. Skip 1st row which is a comment
    # Use the first .counts file in the folder to get the GeneID column as the initial file
    awk -F '\t' '(NR>1)' $(ls -tr "$experiment"/counts/*.counts | head -1) | awk -F '\t' '{print $1}' - > "$experiment"/"$experiment".counts.tsv

    for name in "$experiment"/counts/*.counts
    do
        # Combine the latest test.csv with the count column (last column $NF) of the next file 
        awk -F '\t' '(NR>1)' ${name} | awk -F '\t' '{print $NF}' - | paste -d'\t' "$experiment"/"$experiment".counts.tsv - > temp.tsv
        # Re-write the test.tsv file with the count data from the last file added
        cp temp.tsv "$experiment"/"$experiment".counts.tsv
    done

    rm temp.tsv
    
    cp "$experiment"/"$experiment".counts.tsv .
    
done
