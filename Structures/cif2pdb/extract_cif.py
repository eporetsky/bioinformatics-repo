#https://stackoverflow.com/questions/2018512/reading-tar-file-contents-without-untarring-it-in-python-script
import os
import glob

os.system("mkdir ramdisk/cif/")
os.system("mkdir ramdisk/pdb/")

print("Start extracting cif.gz from tar files")
for tar in glob.glob("*.tar"):
    os.system("tar --wildcards -C ramdisk/cif/ -xvf {0} '*.cif.gz'".format(tar))

print("Start decompressing all cig.gz files")
# Too many files for gunzip to handle by itself
for cif	in glob.glob("ramdisk/cif/*.cif.gz"):
    os.system("gunzip {0}".format(cif)) 
