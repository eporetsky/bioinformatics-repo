import glob, os
import pandas as pd

########################## Conda environment ###########################
# conda install -c conda-forge libgcc-ng=11.2.0
# conda install -c bioconda fastp=0.23.2
# conda install -c bioconda hisat2=2.2.1
# conda install -c bioconda sambamba=0.8.2
# conda install -c bioconda parallel-fastq-dump
# conda install -c bioconda subread
# conda config --add channels conda-forge
# conda install -c conda-forge genozip

########################## Globals ###########################
# genome_index_name: I renamed the hisat2, gtf and genozip genome files to the genome_index_name+suffix
# All these files are contained in the global fil

########################## Running workflow ###########################
# For RNAseq workflow testing change the fastp --reads_to_process to
# a small number (10000) for faster and simpler workflow testing (only applies to fastp, not download step)

# To prolong the life of the hard-drive I recommend to mount ramdisk on the "folder" variable (SRA project name)
# mkdir PRJNA123456
# sudo mount -t tmpfs -o size=80000m tmpfs PRJNA123456/

# To run this snakemake workflow use:
# snakemake -p -j1 --snakefile sra2counts_ena.smk --config sra=PRJNA123456 idx=genome_index_name

def get_reports(acc):
    # Based on the EBI ENA tutorial on accessing their data
    # https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html
    # Change Study Accession (PRJNAXXXXXX) to get the 
    # Change "fields=all" to get all columns in the metadata table
    res = "read_run"
    fie = "study_accession,tax_id,scientific_name,instrument_model,library_strategy,read_count,run_alias,sample_alias,fastq_ftp,fastq_md5"

    # Generate the url that will containt the study accession information once opened
    url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession={0}&result={1}&fields={2}".format(acc, res, fie)
    # Open the URL with pandas and load it as a dataframe
    return pd.read_csv(url, sep="\t")

accession_df = get_reports(config['sra'])
accession_df = accession_df.head()
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
    expand("counts/{sample}.counts", sample=SAMPLES),
    expand("genozip/{sample}.genozip", sample=SAMPLES),

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
      shell("axel -a -n 32 -o ramdisk {params.ftp1}")
      shell("axel -a -n 32 -o ramdisk {params.ftp2}")

  # hisat2 runs faster when using non-compressed fastq files
  rule run_fastp:
    input:
      "ramdisk/{sample}_1.fastq.gz",
      "ramdisk/{sample}_2.fastq.gz",
    output:
      temp("ramdisk/{sample}_1.fastq"),
      temp("ramdisk/{sample}_2.fastq"),
      "reports/{sample}.html"
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
      temp("ramdisk/{sample}_1.fastq"),
      temp("ramdisk/{sample}_2.fastq"),
    output:
      "ramdisk/{sample}.bam"
    params:
      idx = genome_index
    priority: 3
    shell:
      """hisat2 -p 32 --max-intronlen 6000 -x {params.idx} -1 {input[0]} -2 {input[1]} | \
      sambamba view -S -f bam -o /dev/stdout /dev/stdin | \
      sambamba sort  --tmpdir="tmpmba" -t 32 -o {output} /dev/stdin
      """

  rule run_featureCounts_genozip:
    input:
        "ramdisk/{sample}.bam",
    output:
        "counts/{sample}.counts",
        "genozip/{sample}.genozip",
    params:
        idx = genome_index,
    priority: 4
    run:
        shell("featureCounts -p -t exon,CDS -T 32 -a {params.idx}.gtf -o {output[0]} {input}")
        shell("genozip -e {params.idx}.ref.genozip -i bam -o {output[1]} {input}")