# Working with Predicted Structure Data

## 1. cif2pdb

A couple of scripts/attempts to convert large batches of Alphafold2 CIF files to PDB files.

## 2. pLM

A slightly modified notebook obtained from the Rost lab that used ProtT5 to generate protein
and residue embeddings, plus some predictions, for large files. Takes about an hour to run
40k proteins with my RTX2070 on CUDA. I ran a similar script on a workstation with 32 CPU
threads and it took about 30 hours to complete.