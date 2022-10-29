# Phylogenomics Related

## 1. Prior to generating the tree install the following:

```
# FAMSA v2.2.2, ClipKit v1.3., IQtree v2.2.0.3
conda install -c bioconda famsa, clipkit, iqtree

# I also use figtree to visualize the tree
# http://tree.bio.ed.ac.uk/software/Figtree/
# https://github.com/rambaut/figtree/
# On linux an mac change to the lib folder and run:
java -jar ./figtree.jar
```

## 2. For a simple tree building procedure I often use this workflow:
```
# 1. Generate the alignment file:
famsa seqs.fa seqs.aln

# 2. Trim the alignment prior to tree building:
clipkit seqs.aln

# 3. Combine ModelFinder, tree search, ultrafast bootstrap and SH-aLRT test:
iqtree2 -T 32 -s seqs.aln.ckipkit --alrt 1000 -B 1000
```

## 3. Additional scripts:

### GeneID2symbol.py: Convert node names in a .treefile
```
Expects "symbol_converter.csv" and "iqtree.treefile" files and output "iqtree.symbols.treefile"
```

### xml_svg_tree_color_changer.py: Change the color of a node text inside an svg file
```
# Need to find if I have a proper working version of this python script
```