# SUPERFAMILY Proteome Annotation
You can find the SUPERFAMILY database of structural and functional annotation here: https://supfam.org/SUPERFAMILY/\
To download the stand-alone SUPERFAMILY scripts you will need to register the licence: https://supfam.org/SUPERFAMILY/register.html\
There are O.K. installation and running instructions in the different tabs here: https://supfam.org/SUPERFAMILY/howto_use_models.html\

The SCOP links in the tutorial are broken but you can use this instead\
```
# Source: http://scop.mrc-lmb.cam.ac.uk/download
wget http://scop.mrc-lmb.cam.ac.uk/files/scop-cla-latest.txt
wget http://scop.mrc-lmb.cam.ac.uk/files/scop-des-latest.txt
# I don't think I ended up using the scop-des for anything though.
```

If you installed the latest version of HMMER (>v3.2) you will get an error message when using hmmpress\
You can use the hmmer-3.1b2.tar. I will include the built hmmpress (v3.1b2) in this repository just in case.\

I have replaced all $file with $ARGV[0] because it can't handle "complicated" file names otherwise.\
This will keep the full fasta file name, but I didn't feel like learning perl regex for this.\

I like to run the pipeline on ramdisk as well. Reduce read-write of large intermediate files.\
I have included a short bash script that automatically unzips and runs the SUPERFAMILY pipeline.\
```
mkdir ramdisk
sudo mount -t tmpfs -o size=80000m tmpfs ramdisk/
cd ramdisk
for d in *.gz ; do ( gunzip --keep --force $d && ./superfamily.pl ${d%%.gz}); done
```

Additional short python script to get the name of the SUPERFAMILY protein family in the result file\
```
import pandas as pd
models = pd.read_csv("model.tab", sep="\t", header=None, dtype=str)
ass = pd.read_csv("results.fa.ass", sep="\t", header=None)
merged = pd.merge(ass, models, left_on=1, right_on=0)
merged = merged[["0_x", "1_x", "3_x", "4_y"]]
merged.columns = ["gene_ID", "SSF", "eval", "description"]
merged.to_csv("results.ssf", sep="\t")
```
