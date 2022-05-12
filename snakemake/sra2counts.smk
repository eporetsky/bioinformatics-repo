import glob, os
import pandas as pd

########################## Conda environment ###########################
# conda config --add channels conda-forge
# conda install genozip
# conda install -c conda-forge libgcc-ng=11.2.0
# conda install -c bioconda fastp=0.23.2
# conda install -c bioconda hisat2=2.2.1
# conda install -c bioconda sambamba=0.8.2
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
# snakemake -F -p -j 1 --snakefile fastq2genozip.smk --config sra=PRJNA123456 idx=genome_index_name --force


try:
  folder = config['sra']
  os.system("wget 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term="+folder+"'"+' -O - | tee '+folder+'.csv')
except:
  folder = ''

SAMPLES = list(pd.read_csv(folder+'.csv')["Run"])

fastq_type = pd.read_csv(folder+'.csv')["LibraryLayout"][0]

# If reads_to_process not defined, set to 0 to process all reads
try: rtp = config['rtp'] 
except: rtp=0

try: genome_index = "global/"+config['idx']
except: genome_index = glob.glob("global/*.1.ht2")[0].replace(".1.ht2", "")

# If a the folder was not previously created (and ideally mounted as a ramdisk)
try:
  os.mkdir(folder)
except:
  None

try:
  os.mkdir(folder+"/fastp")
  os.mkdir(folder+"/mapped")
  os.mkdir(folder+"/reports") 
  os.mkdir("counts")
  os.mkdir("genozip")
except:
  None

########################## Rules ###########################

rule all:
  input:
    # The snakemake workflow should include these files
    expand("counts/{sample}.counts", sample=SAMPLES),
    expand("genozip/{sample}.genozip", sample=SAMPLES),

# Run fastp trimming on the fastq file
if fastq_type == "PAIRED":
  rule run_fastq_dump:
    output:
      temp(folder+"/{sample}_1.fastq.gz"),
      temp(folder+"/{sample}_2.fastq.gz"),
    params:
      fld = folder,
    priority: 1
    shell:
      """parallel-fastq-dump --sra-id {wildcards.sample} --threads 32 --outdir {params.fld} --tmpdir {params.fld} --split-files --gzip"""

  # hisat2 runs faster when using non-compressed fastq files
  rule run_fastp:
    input:
      folder+"/{sample}_1.fastq.gz",
      folder+"/{sample}_2.fastq.gz",
    output:
      temp(folder+"/fastp/{sample}_1.fastq"),
      temp(folder+"/fastp/{sample}_2.fastq"),
      folder+"/reports/{sample}.html"
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
      folder+"/fastp/{sample}_1.fastq",
      folder+"/fastp/{sample}_2.fastq",
    output:
      temp(folder+"/mapped/{sample}.bam")
    params:
      idx = genome_index
    priority: 3
    shell:
      """hisat2 -p 32 --max-intronlen 6000 -x {params.idx} -1 {input[0]} -2 {input[1]} | \
      sambamba view -S -f bam -o /dev/stdout /dev/stdin | \
      sambamba sort  --tmpdir="tmpmba" -t 32 -o {output} /dev/stdin
      """

  rule run_featureCounts:
    input:
        folder+"/mapped/{sample}.bam",
    output:
        "counts/{sample}.counts",
    params:
        idx = genome_index,
    priority: 4
    shell:
        "featureCounts -p -t exon,CDS -T 32 -a {params.idx}.gtf -o {output} {input}"

  rule run_genozip:
    input:
        folder+"/mapped/{sample}.bam",
    output:
        "genozip/{sample}.genozip",
    params:
        idx = genome_index,
    priority: 5
    shell:
        "genozip -e {params.idx}.genozip -i bam -o {output} {input}"