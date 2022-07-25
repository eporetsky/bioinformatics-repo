# Renaming protein sequence fasta files using a csv mapping file
Because OrthoFinder (and other tools) use the fasta name for species/genotype name\
I have added the script below to take the full portein fasta file name and convert it\
to the preferred name. Keep the mapping file to remember which version was used.
```
# https://stackoverflow.com/questions/50263422/renaming-files-from-csv-mapping-using-bash-mangles-files
sed 's/"//g' name_mapping.csv | while IFS=, read orig new; do mv "$orig" "$new"; done

# csv should containt one row for each file: "file_name.fa, species_name.fa"
```

# Removing isoform numbers from longest/primary/canonical protein fasta files

```
Use python script found in: https://github.com/eporetsky/bioinformatics-assortment/tree/master/Genomes/Scripts/Canonical

For a single file:
python clean_fasta_isoforms.py Zmays_493_RefGen_V4.protein_primaryTranscriptOnly.primary.primary.fa

For all .fa files in folder (will need to edit python script for all .fasta files)
python clean_fasta_isoforms.py

The folder contains additional scripts that help extract the canonical/primary transcripts for each fasta file\
with some database- and genome-specific implemntations.
```

# Run a simple Diamond Top Blast Hit

```
# Subject will be in the first column and include all genes, even with no hits
diamond makedb --in query.fa --db query
diamond blastp -d query -q subject.fa -o results.tsv --max-target-seqs 1 --unal 1 --quiet --header
```

# Run OrthoFinder

```
orthofinder -t 32 -a 32 -f seqs/
```

# Run Broccoli

```
python broccoli.py -dir seqs -threads 32
```
