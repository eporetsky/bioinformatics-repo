#https://stackoverflow.com/questions/2018512/reading-tar-file-contents-without-untarring-it-in-python-script
import os
import glob
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from Bio.PDB import MMCIFParser, PDBIO

os.system("mkdir ramdisk/cif/")
os.system("mkdir ramdisk/pdb/")

print("Start extracting cif.gz from tar files")
for tar in glob.glob("*.tar"):
    os.system("tar --wildcards -C ramdisk/cif/ -xvf {0} '*.cif.gz'".format(tar))

print("Start decompressing all cig.gz files")
# Too many files for gunzip to handle by itself
for cif	in glob.glob("ramdisk/cif/*.cif.gz"):
    os.system("gunzip {0}".format(cif)) 

print("Start convert cif to pdb files")
# Convert all saved cif.gz to pdb format
# This step is very slow on a single node
# Re-implemented as a parallel script
def cif2pdb(cif_path):
    parser = MMCIFParser()
    # The parser can handle both compressed and not so can just replace text for either case
    struct_id = cif_path.replace("ramdisk/cif/", "").replace(".cif", "")
    # Read the CIF file and give it the base filename
    struct = parser.get_structure(struct_id, cif_path)
    # Convert the CIF structure into a PDB structure and save
    pdbio = PDBIO()
    pdbio.set_structure(struct)
    pdbio.save("ramdisk/pdb/"+struct_id+".pdb")

# Run the cif2pdb conversion in parallel
pool = ThreadPool(32)
pool.map(cif2pdb, glob.glob("ramdisk/cif/*.cif"))
