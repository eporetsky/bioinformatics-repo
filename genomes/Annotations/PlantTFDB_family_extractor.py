import pandas as pd
from Bio import SeqIO

# Which TF family sequences to extract
tf_family = "WRKY"

# Load the TF annotation output from PlantTFDB
# http://planttfdb.gao-lab.org/prediction.php
df = pd.read_csv("TFs_Athaliana.tsv", sep="\t", header=None)
df = df[df[1].isin([tf_family])]
gene_list = df[0].to_list()

# To split the gene ID column and remove isoform IDs
# df[0].str.split(".",expand=True)[0].to_list()

with open("At"+tf_family+".fasta","w") as f:
    for gene in gene_list:
        # The sequences I've used were downloaded from Phytozome
        for seq_record in SeqIO.parse("Athaliana_167_TAIR10.protein_primaryTranscriptOnly.fa", "fasta"):
            if gene == seq_record.id:
                # Write the GeneID to the file. I remove isoform IDs. Adjust accordingly
                f.write(">" + str(seq_record.id).split(".")[0] + "\n") # for Arabidopsis
                ### f.write(">" + str(seq_record.id).split("_")[0] + "\n") # for Maize
                ### f.write(">" + ".".join(str(seq_record.id).split(".")[0:2]) + "\n") # for Sorghum
                
                # Write the sequence to the fasta file. Check if ends with * to remove it
                if str(seq_record.seq)[-1] == "*":
                    f.write(str(seq_record.seq)[:-1] + "\n") 
                else:
                    f.write(str(seq_record.seq) + "\n")