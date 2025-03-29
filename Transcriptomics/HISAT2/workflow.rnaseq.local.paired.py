# sudo mount -t tmpfs -o size=250000m tmpfs ramdisk/

# sudo mount -t tmpfs -o size=80000m tmpfs ramdisk/
import sys
import os
import glob
import pandas as pd

experiment = sys.argv[1]
genome = sys.argv[2]

os.makedirs(os.path.join(experiment, "crams"), exist_ok=True)
os.makedirs(os.path.join(experiment, "counts"), exist_ok=True)
os.makedirs(os.path.join(experiment, "reports"), exist_ok=True)

fl_df = []
for fl in glob.glob(f"{experiment}/fastq/*_1.fq.gz"):
    fl = os.path.basename(fl)
    name = fl.replace("_1.fq.gz", "")
    fl_df.append([experiment, name, genome])

fl_df = pd.DataFrame(fl_df)
fl_df.to_csv("process.tsv", sep="\t", index=None, header=None)
os.system("parallel -j1 --colsep '\\t' -a process.tsv bash workflow.rnaseq.local.paired.sh {1} {2} {3}")