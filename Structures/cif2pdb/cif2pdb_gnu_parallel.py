#https://stackoverflow.com/questions/2018512/reading-tar-file-contents-without-untarring-it-in-python-script

import sys
from Bio.PDB import MMCIFParser, PDBIO

cif_path = sys.argv[1]
parser = MMCIFParser()
# The parser can handle both compressed and not so can just replace text for either case
struct_id = cif_path.replace("cif/", "").replace(".cif", "")
# Read the CIF file and give it the base filename
struct = parser.get_structure(struct_id, "cif/"+cif_path)
# Convert the CIF structure into a PDB structure and save
pdbio = PDBIO()
pdbio.set_structure(struct)
pdbio.save("pdb/"+struct_id+".pdb")
