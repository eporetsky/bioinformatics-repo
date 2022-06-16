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

########################## Globals ###########################
# genome_index_name: I renamed the hisat2, gtf and genozip genome files to the genome_index_name+suffix
# All these files are contained in the global fil

########################## Running workflow ###########################
# For RNAseq workflow testing change the fastp --reads_to_process to
# a small number (10000) for faster and simpler workflow testing (only applies to fastp, not download step)

# To prolong the life of the hard-drive I recommend to mount ramdisk on the "folder" variable (SRA project name)
# mkdir PRJNA123456
# sudo mount -t tmpfs -o size=95000m tmpfs PRJNA123456/

# To run this snakemake workflow use:
# snakemake -F -p -j1 --keep-going --snakefile sra2counts_ena.smk --config sra=PRJNA123456 idx=Zmays_493_APGv4_Phytozome

fq_list = []
fastq_type = "SINGLE"
for fl in os.listdir():
    if "fq.gz" in fl or "fastq.gz" in fl :
        fq_list.append(fl)
        if "_2" in fl:
            fastq_type = "PAIRED"
            
if fastq_type=="PAIRED":
    SAMPLES = [fq.split("_")[0] for fq in fq_list] 
else:
    SAMPLES = [fq.split(".")[0] for fq in fq_list]

# If reads_to_process not defined, set to 0 to process all reads
try: rtp = config['rtp'] 
except: rtp=0

try: genome_index = "global/"+config['idx']
except: genome_index = glob.glob("global/*.1.ht2")[0].replace(".1.ht2", "")

########################## Rules ###########################

rule all:
  input:
    # The snakemake workflow should include these files
    expand("counts/{sample}.counts.gz", sample=SAMPLES),
    expand("genozip/{sample}.genozip", sample=SAMPLES),

# Run fastp trimming on the fastq file
if fastq_type == "SINGLE":
  # hisat2 runs faster when using non-compressed fastq files
  rule run_fastp:
    input:
      "ramdisk/{sample}.fastq.gz",
    output:
      temp("ramdisk/{sample}.fastq"),
      "reports/{sample}.html"
    params:
      rtp = rtp,
    priority: 1
    run:
      shell("rm -f -r ramdisk/*")
      shell("""fastp\
        --reads_to_process {params.rtp}\
        --in1 {input} --out1 {output[0]}\
        --html {output[1]}""")
  
  rule run_hisat2:
    input:
      temp("ramdisk/{sample}.fastq"),
    output:
      "reports/hisat2_{sample}.txt",
      "ramdisk/{sample}.bam"
    params:
      idx = genome_index
    priority: 2
    shell:
      """hisat2 -p 32 --max-intronlen 6000 -x {params.idx} -U {input} --summary-file {output[0]} | \
      sambamba view -S -f bam -F "not unmapped" -o /dev/stdout /dev/stdin | \
      sambamba sort  --tmpdir="tmpmba" -t 32 -o {output[1]} /dev/stdin
      """

  rule run_featureCounts_genozip:
    input:
        temp("ramdisk/{sample}.bam"),
    output:
        "counts/{sample}.counts.gz",
        "genozip/{sample}.genozip",
    params:
        idx = genome_index,
        smpl = lambda wc: wc.get("sample"),
    priority: 3
    run:
        shell("featureCounts -t exon,CDS -T 32 -a {params.idx}.gtf -o counts/{params.smpl}.counts {input}")
        shell("gzip counts/{params.smpl}.counts")
        shell("genozip -e {params.idx}.ref.genozip -i bam -o {output[1]} {input}")

# Run fastp trimming on the fastq file
if fastq_type == "PAIRED":

  # hisat2 runs faster when using non-compressed fastq files
  rule run_fastp:
    input:
      "{sample}_1.fastq.gz",
      "{sample}_2.fastq.gz",
    output:
      temp("ramdisk/{sample}_1.fastq"),
      temp("ramdisk/{sample}_2.fastq"),
      "reports/fastp_{sample}.html"
    params:
      rtp = rtp,
    priority: 1

    run:
      shell("rm -f -r ramdisk/*")
      shell("""fastp\
        --reads_to_process {params.rtp}\
        --in1 {input[0]} --out1 {output[0]}\
        --in2 {input[1]} --out2 {output[1]}\
        --html {output[2]}""")
  
  rule run_hisat2:
    input:
      temp("ramdisk/{sample}_1.fastq"),
      temp("ramdisk/{sample}_2.fastq"),
    output:
      "ramdisk/{sample}.bam",
      "reports/hisat2_{sample}.txt",
    params:
      idx = genome_index
    priority: 2
    shell:
      """hisat2 -p 32 --max-intronlen 6000 -x {params.idx} -1 {input[0]} -2 {input[1]} --summary-file {output[1]} | \
      sambamba view -S -f bam -F "not unmapped" -o /dev/stdout /dev/stdin | \
      sambamba sort  --tmpdir="tmpmba" -t 32 -o {output[0]} /dev/stdin
      """

  rule run_featureCounts_genozip:
    input:
        temp("ramdisk/{sample}.bam"),
    output:
        "counts/{sample}.counts.gz",
        "genozip/{sample}.genozip",
    params:
        idx = genome_index,
        smpl = lambda wc: wc.get("sample"),
    priority: 3
    run:
        shell("featureCounts -p -t exon,CDS -T 32 -a {params.idx}.gtf -o counts/{params.smpl}.counts {input}")
        shell("gzip counts/{params.smpl}.counts")
        shell("genozip -e {params.idx}.ref.genozip -i bam -o {output[1]} {input}")