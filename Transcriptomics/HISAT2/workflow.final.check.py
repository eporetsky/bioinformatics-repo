# sudo mount -t tmpfs -o size=180000m tmpfs ramdisk/

import glob, os

for acc in glob.glob("*/"):
    acc = acc.replace("/", "")
    fl_list = list(glob.glob(f"{acc}/*/"))
    if f"{acc}/crams/" not in fl_list:
        continue
    crams_set = set([os.path.basename(fl).replace(".cram", "") for fl in glob.glob(f"{acc}/crams/*")])
    fastp_set = set([os.path.basename(fl).replace(".html", "").replace("fastp_", "") for fl in glob.glob(f"{acc}/reports/*.html")])
    hisat_set = set([os.path.basename(fl).replace(".txt", "").replace("hisat2_", "") for fl in glob.glob(f"{acc}/reports/*.txt")])
    count_list = set([os.path.basename(fl).replace(".counts", "") for fl in glob.glob(f"{acc}/counts/*.counts")])
    
    miss_fastp = crams_set.difference(fastp_set)
    miss_hisat = crams_set.difference(hisat_set)
    miss_count = crams_set.difference(count_list)

    all_miss = miss_fastp.union(miss_hisat).union(miss_count)
    
    for missing in list(all_miss):
        os.system(f"rm {acc}/crams/{missing}.cram")
