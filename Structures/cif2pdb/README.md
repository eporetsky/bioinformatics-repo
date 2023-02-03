# Working with Predicted Structure Data

## Step 1

Have all the tar files in a folder and make the following folders:
```
mkdir ramdisk
# if you want to mount the ramdisk folder as a ramdisk then use:
sudo mount -t tmpfs -o size=95000m tmpfs ramdisk/

# also make these folders:
mkdir ramdisk/cif
mkdir ramdisk/pdb

```

Next, extract all the CIF files from all the TAR files in the current folder:
```
python extract_cif
```

Then convert all the CIF files to PDB (might take a while, specially on few threads)
```
cd ramdisk
ls cif/ | parallel "python cif2pdb_gnu_parallel.py  {}"
```

Finally, move all the CIF and PDB file to separate TAR files
```
cd pdb
tar -czf ../name.pdb.tar.gz .`
cd ../cif
tar -czf ../name.cif.tar.gz .`
```

I think this should work to decompress the TAR files inside a folder
```
tar cvf - name.cif.tar.gz
tar cvf - name.pdb.tar.gz
```

# Older scripts

## 1. cif2pdb.py

This python script is meant to parse AlphaFold2 non-model organism tar files, 
extract all the cif.gz files and convert them to pdb format for downstream analyses.
It works but it rather slow.

## 2. cif2pdb_parallel_ramdisk.py

Does the same as above but I tried to use the multiprocessing python library to
speed it up. Still pretty slow and I wouldn't recommend using for more than 10k CIF files.

