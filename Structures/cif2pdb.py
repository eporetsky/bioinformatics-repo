#https://stackoverflow.com/questions/2018512/reading-tar-file-contents-without-untarring-it-in-python-script
import os
import glob
from Bio.PDB import MMCIFParser, PDBIO

os.system("mkdir cif/")
os.system("mkdir pdb/")

for tar in glob.glob("*.tar"):
    os.system("tar -C cif/ -xvf {0} '*.cif.gz'".format(tar))

# Too many files for gunzip to handle by itself
for cif	in glob.glob("cif/*.cif.gz"):
    os.system("gunzip {0}".format(cif)) 

# Convert all saved cif.gz to pdb format
for cif_path in glob.glob("cif/*.cif"):    
    parser = MMCIFParser()
    # The parser can handle both compressed and not so can just replace text for either case
    struct_id = cif_path.replace("cif/", "").replace(".cif", "")
    # Read the CIF file and give it the base filename
    struct = parser.get_structure(struct_id, cif_path)
    # Convert the CIF structure into a PDB structure and save
    pdbio = PDBIO()
    pdbio.set_structure(struct)
    pdbio.save("pdb/"+struct_id+".pdb")