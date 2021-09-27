import pandas as pd

# csv file with a GeneID and gene symbol header
# https://stackoverflow.com/questions/17140886/how-to-search-and-replace-text-in-a-file
conversion_dict = pd.read_csv("symbol_converter.csv", index_col="GeneID").to_dict()["Symbol"]

# Read in the file
with open("iqtree.treefile", 'r') as file :
  treefile = file.read()

# Replace the target string
for gene in conversion_dict.keys():
    treefile = treefile.replace(gene, conversion_dict[gene])

# Write the file out again
with open('iqtree.symbols.treefile', 'w') as file:
  file.write(treefile)