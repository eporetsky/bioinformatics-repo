import os
import glob
import pandas as pd

def extract_gene_ids(fasta_file):
    gene_ids = []
    
    with open(fasta_file, 'r') as file:
        for line in file:
            # Check if the line starts with '>'
            if line.startswith('>'):
                # Extract the gene ID from the header line
                gene_id = line[1:].strip().split()[0]  # Get the first part after '>'
                gene_ids.append(gene_id)
    
    return gene_ids

genome_dict = {}
for fl in glob.glob("clean/*"):
    genome_id = os.path.basename(fl).replace(".fa", "")
    genome_dict[genome_id] = extract_gene_ids(fl)

#annots = pd.DataFrame()
tmp_list = []
for fl in glob.glob("output/*"):
    tmp = pd.read_csv(fl, sep="\t", header=None)
    tmp_list.append(tmp)
annots = pd.concat(tmp_list)
annots.columns = ['gene', 'transcript_id', 'length', 'database', 'domain_id', 'domain_name', 
           'start', 'end', 'evalue', 'status', 'date', 'entry_type', 'entry_name', "go", "reactome_metacyc"]
annots = annots.drop("reactome_metacyc", axis=1)

mappings = pd.read_csv("sequence_mapping.tsv", sep="\t")
mappings.columns = ["geneID", "genome_id", "gene"]

os.makedirs("annots", exist_ok=True)
for fl in glob.glob("clean/*"):
    print(fl)
    genome_id = os.path.basename(fl).replace(".fa", "")
    #tmp_df = mappings[mappings["geneID"].isin(genome_dict[genome_id])]
    tmp_df = mappings[mappings["genome_id"]==genome_id]
    tmp_df = tmp_df.drop("genome_id", axis=1)
    annot_df = tmp_df.merge(annots, on="gene", how="left")
    annot_df = annot_df.dropna()
    annot_df = annot_df.drop("gene", axis=1)
    annot_df.to_csv(f"annots/{genome_id}.ips.tsv.gz", sep="\t", index=False)