# Filtering non-biallelic sites from a hapmap file

rMVP requries the files to only contain biallelic sites. The script below keeps only the common single nucleotide SNPs.

```
awk '{ if ($2 == "A/T" || $2 == "A/G" || $2 == "A/C" || \
           $2 == "T/A" || $2 == "T/G" || $2 == "T/C" || \
           $2 == "G/A" || $2 == "G/T" || $2 == "G/C" || \
           $2 == "C/A" || $2 == "C/G" || $2 == "C/T" ||$2 == "alleles") { print } }' snps.hmp.txt | snps.biallelic.hmp.txt
```