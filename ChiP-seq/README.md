# ChiP-seq Analysis Workflow
Just started working on this and I will make updates as soon as I can.\

Currently working with Arabidopsis data.\ 
Genome was downloaded from Phytozome Araport11 (Athaliana_447_Araport11).\


## Clean genome fasta file
```
# The fasta headers contain ';' characters which is not compatible with Genozip. Bash script to clean is based on link below
# https://yedomon.netlify.app/posts/2021-01-12-awk-how-to-remove-the-rest-of-a-fasta-header-name-after-a-specific-character/
cat Athaliana_447.fa | awk -F' ' '{print $1}' > Athaliana_447.clean.fa
```


## Prepare indexes for processing and compressing files
```
# Build hisat2 genome index for aligning
hisat2-build -p 32 Athaliana_447_TAIR10.fa Athaliana_447

# Create a gtf file (not sure this will be needed but just in case)
# gffread Athaliana_447_Araport11.gene.gff3 -T -o Athaliana_447.gtf

# Not critical but I will eventually compress outputs using genozip
genozip --make-reference Athaliana_447.clean.fa
```
