# sudo mount -t tmpfs -o size=180000m tmpfs ramdisk/

import glob, os
import pandas as pd

experiment_list = []
for experiment in os.listdir("."):
    if "ramdisk" in experiment:
        continue
    # check whether the current object is a folder or not
    if os.path.isdir(os.path.join(os.path.abspath("."), experiment)): 
        experiment_list.append(experiment)
    if ".gtf" in experiment:
        index_name = experiment.split(".gtf")[0]

for experiment in experiment_list:
    os.makedirs(os.path.join(experiment, "crams"), exist_ok=True)
    os.makedirs(os.path.join(experiment, "counts"), exist_ok=True)
    os.makedirs(os.path.join(experiment, "reports"), exist_ok=True)

    url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession={0}&result=read_run&fields=fastq_ftp".format(experiment)
    # Open the URL with pandas and load it as a dataframe
    try:
        accession_df = pd.read_csv(url, sep="\t")
    except:
        print("Did not work:", experiment)
        continue
    
    accession_df["num_fq"] = accession_df["fastq_ftp"].apply(lambda x: len(x.split(";")))

    done_list = []
    for fl in glob.glob(os.path.join(experiment, "crams", "*")):
        done_list.append(os.path.splitext(os.path.basename(fl))[0])

    paired_df = []
    single_df = []
    for key, val in accession_df[["run_accession", "fastq_ftp"]].values.tolist():
        splits = val.split(";")
        # Sometimes there are paired + single fastq files, keep only paired
        if len(splits)==3:
            splits = splits[1:]
        if key not in done_list:
            if len(splits)==2:
            	paired_df.append([experiment, key, splits[0], splits[1], index_name])
            if len(splits)==1:
            	single_df.append([experiment, key, splits[0], index_name])
            	
    print(f"Paired: {len(paired_df)} and Single: {len(single_df)}")
    
    if len(paired_df)>0:
    	paired_df = pd.DataFrame(paired_df)
    	paired_df.to_csv("process.tsv", sep="\t", index=None, header=None)
    	print("Starting to run {}. Done: {} Remaining: {}".format(experiment, len(done_list), len(paired_df)+len(single_df)))
    	os.system("parallel -j2 --colsep '\\t' -a process.tsv bash workflow.rnaseq.paired.double.sh {1} {2} {3} {4} {5}")
    
    if len(single_df)>0:
    	single_df = pd.DataFrame(single_df)
    	single_df.to_csv("process.tsv", sep="\t", index=None, header=None)
    	print("Starting to run {}. Done: {} Remaining: {}".format(experiment, len(done_list), len(paired_df)+len(single_df)))
    	os.system("parallel -j2 --colsep '\\t' -a process.tsv bash workflow.rnaseq.single.double.sh {1} {2} {3} {4}")
    	#os.system("parallel -j3 --colsep '\\t' -a process.tsv echo {1} {2} {3} {4}")
