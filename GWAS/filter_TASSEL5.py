import pandas as pd

# Filter TASSEL5 statistics output
df = pd.read_csv("gwas_results.txt", sep="\t")

# Keep columns needed for Manhattan plot
df = df[["Trait", "Chr", "Pos", "p"]]

# Filter by p-value and ascending sort column
# df[df["p"]<0.0001].sort_values(by="p", ascending=True)

# Save file as tsv
df.to_csv("gwas_results.mini.tsv", sep="\t")
