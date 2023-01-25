# Working with Predicted Structure Data

## 1. cif2pdb.py

This python script is meant to parse AlphaFold2 non-model organism tar files, 
extract all the cif.gz files and convert them to pdb format for downstream analyses.
It works but it rather slow.

## 2. cif2pdb_parallel_ramdisk.py

Does the same as above but I tried to use the multiprocessing python library to
speed it up. Still pretty slow and I wouldn't recommend using for more than 10k CIF files.

## 3. cif2pdb_gnu_parallel.py

Uses the linux GNU parallel package (`conda install -c conda-forge parallel`).
Makde a pdb directory and use the output cif/ folder from cif2pdb_parallel_ramdisk.py
`ls cif/ | parallel "python cif2pdb_gnu_parallel.py  {}"`