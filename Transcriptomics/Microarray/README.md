# Converting Genbank IDs from microarry metadata to current gene IDs

## Steps 
* Obtain GB ids
* Download GB sequences 
* Blast sequences against all genes
* Remove duplicates

## To generate a genozip reference genome file
```
Use: https://github.com/LeeBergstrand/Genbank-Downloaders
# Download GetNucleotide.py and SeqExtract.py
# Create a txt file that contains a single column with all GB ids
# Run in an empty folder because it generates a file for each sequence
python GetNucleotide.py updated_gb.txt your@email.edu
```

## Generate a blast db for the mRNA or cDNA fasta file for your genome
```
gunzip -c -k Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gene.fa.gz | makeblastdb -title "B73v4" -dbtype nucl -in --out B73v4
```

## Run blast on all downloaded sequences (if possible run in parallel) 
```
# Find all .fna files in subdirectory
find -mindepth 2 -maxdepth 2 -name *.fna | xargs -n1 -P30  -I{} sh -c 'blastn -db B73v4 -query {} -outfmt "6 qseqid sseqid evalue bitscore" >> blast_results_parallel.txt'
```

## Sort by e-value or bit score and remove duplicates
```
I used excel. Will use R/py if I run larger files
```
