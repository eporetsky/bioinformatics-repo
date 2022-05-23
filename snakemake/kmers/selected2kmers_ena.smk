import glob, os
import pandas as pd

########################## Conda environment ###########################
# conda config --add channels conda-forge
# conda install genozip
# conda install -c conda-forge libgcc-ng=11.2.0
# conda install -c bioconda fastp=0.23.2
# conda install -c bioconda hisat2=2.2.1
# conda install -c bioconda sambamba=0.8.2
# conda install -c bioconda parallel-fastq-dump
# conda install -c conda-forge genozip
# sudo apt-get install axel

########################## Globals ###########################
# genome_index_name: I renamed the hisat2, gtf and genozip genome files to the genome_index_name+suffix
# All these files are contained in the global fil

########################## Running workflow ###########################
# For RNAseq workflow testing change the fastp --reads_to_process to
# a small number (10000) for faster and simpler workflow testing (only applies to fastp, not download step)

# To prolong the life of the hard-drive I recommend to mount ramdisk on the "folder" variable (SRA project name)
# mkdir ramdisk
# sudo mount -t tmpfs -o size=95000m tmpfs ramdisk/

# To run this snakemake workflow use:
# snakemake -F -p -j1 --keep-going --snakefile selected2kmers_ena.smk --config tbl=PRJNA123456_selected.csv

accession_df = pd.read_csv(config['tbl'])
fastq_type = "PAIRED" if len(accession_df["fastq_ftp"].values.tolist()[0].split(";")) == 2 else "SINGLE"
accession_dict = {}
SAMPLES = []
for key, val in accession_df[["run_accession", "fastq_ftp"]].values.tolist():
    accession_dict[key] = val.split(";")
    SAMPLES.append(key)

# If reads_to_process not defined, set to 0 to process all reads
try: rtp = config['rtp'] 
except: rtp=0

try: genome_index = "global/"+config['idx']
except: genome_index = glob.glob("global/*.1.ht2")[0].replace(".1.ht2", "")

########################## Rules ###########################

rule all:
  input:
    # The snakemake workflow should include these files
    expand("kmc/{sample}.kmc_suf.gz", sample=SAMPLES),
    expand("kmc/{sample}.kmc_pre.gz", sample=SAMPLES),

# Run fastp trimming on the fastq file
if fastq_type == "SINGLE":
  rule run_fastq:
    output:
      "ramdisk/{sample}.fastq.gz",
    params:
      ftp = lambda wc: accession_dict[wc.get("sample")][0],
    priority: 1
    run:
      shell("rm -f -r ramdisk/*")
      shell("axel --quiet -n 32 -o ramdisk {params.ftp}")

  # hisat2 runs faster when using non-compressed fastq files
  rule run_fastp:
    input:
      "ramdisk/{sample}.fastq.gz",
    output:
      "ramdisk/{sample}.fastq",
      "reports/{sample}.html"
    params:
      rtp = rtp,
    priority: 2
    shell:
      """fastp\
      --reads_to_process {params.rtp}\
      --in1 {input} --out1 {output[0]}\
      --html {output[1]}"""

  rule run_kmc:
    input:
      "ramdisk/{sample}.fastq"
    output:
      "ramdisk/{sample}.kmc_suf",
      "ramdisk/{sample}.kmc_pre"
    priority: 3
    params:
      smpl = lambda wc: wc.get("sample"),
    run:
      shell("kmc -k25 -m100 -ci2 -r {input} ramdisk/{params.smpl} /tmp/")
  
  rule run_gzip:
    input:
      "ramdisk/{sample}.kmc_suf",
      "ramdisk/{sample}.kmc_pre"
    output:
      "kmc/{sample}.kmc_suf.gz",
      "kmc/{sample}.kmc_pre.gz"
    priority: 4
    run:
      shell("gzip {input[0]}")
      shell("mv {input[0]}.gz kmc/")
      shell("gzip {input[1]}")
      shell("mv {input[1]}.gz kmc/")


# Run fastp trimming on the fastq file
if fastq_type == "PAIRED":
  rule run_fastq:
    output:
      temp("ramdisk/{sample}_1.fastq.gz"),
      temp("ramdisk/{sample}_2.fastq.gz"),
    params:
      ftp1 = lambda wc: accession_dict[wc.get("sample")][0],
      ftp2 = lambda wc: accession_dict[wc.get("sample")][1],
    priority: 1
    run:
      shell("rm -f -r ramdisk/*")
      shell("axel --quiet -n 32 -o ramdisk {params.ftp1}")
      shell("axel --quiet -n 32 -o ramdisk {params.ftp2}")

  # hisat2 runs faster when using non-compressed fastq files
  rule run_fastp:
    input:
      "ramdisk/{sample}_1.fastq.gz",
      "ramdisk/{sample}_2.fastq.gz",
    output:
      temp("ramdisk/{sample}_1.fastq"),
      temp("ramdisk/{sample}_2.fastq"),
      "reports/fastp_{sample}.html"
    params:
      rtp = rtp,
    priority: 2

    shell:
      """fastp\
      --reads_to_process {params.rtp}\
      --in1 {input[0]} --out1 {output[0]}\
      --in2 {input[1]} --out2 {output[1]}\
      --html {output[2]}"""
  
  rule run_hisat2:
    input:
      "ramdisk/{sample}_1.fastq",
      "ramdisk/{sample}_2.fastq",
    output:
      temp("ramdisk/{sample}.bam"),
      "reports/hisat2_{sample}.txt",
    params:
      idx = genome_index
    priority: 3
    shell:
      """hisat2 -p 32 --max-intronlen 6000 -x {params.idx} -1 {input[0]} -2 {input[1]} --summary-file {output[1]} | \
      sambamba view -S -f bam -F "not unmapped" -o /dev/stdout /dev/stdin | \
      sambamba sort  --tmpdir="tmpmba" -t 32 -o {output[0]} /dev/stdin
      """

  rule run_featureCounts_genozip:
    input:
        "ramdisk/{sample}.bam",
    output:
        "counts/{sample}.counts.gz",
        "genozip/{sample}.genozip",
    params:
        idx = genome_index,
        smpl = lambda wc: wc.get("sample"),
    priority: 4
    run:
        shell("featureCounts -p -t exon,CDS -T 32 -a {params.idx}.gtf -o counts/{params.smpl}.counts {input}")
        shell("gzip counts/{params.smpl}.counts")
        shell("genozip -e {params.idx}.ref.genozip -i bam -o {output[1]} {input}")