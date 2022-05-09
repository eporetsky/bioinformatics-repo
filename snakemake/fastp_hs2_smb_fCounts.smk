import glob, os

# For RNAseq workflow testing change the fastp --reads_to_process to
# a small number (10000) for faster and simpler workflow testing

# To run this snakemake workflow use:
# single-end: --config sp="single"
# paired-end: --config sp="paired"
# snakemake --snakefile fastp_hs2_smb_fCounts.smk --config sp="single" idx="Zmays_493_APGv4_Phytozome" rtp=10000  --cores 32 --force


########################## Globals ###########################

# https://stackoverflow.com/questions/45971333/use-of-argparse-in-snakemake-script

# conda config --add channels conda-forge
# conda install genozip
# conda install -c conda-forge libgcc-ng=11.2.0
# conda install -c bioconda fastp=0.23.2
# conda install -c bioconda hisat2=2.2.1
# conda install -c bioconda sambamba=0.8.2

fastq_type = config['sp']

if len(glob.glob("*.fastq.gz")) > 0:
  fastq_suffix = "fastq.gz"
elif len(glob.glob("*.fq.gz") > 0 ):
  fastq_suffix = "fq.gz"
else:
  print("Uncertain about fastq suffix")

if fastq_type == "paired":
  SAMPLES = [fl.replace("_1.fastq.gz","") for fl in glob.glob("*_1.fastq.gz")]
elif fastq_type == "single":
  SAMPLES = [fl.replace(".fastq.gz","") for fl in glob.glob("*.fastq.gz")]

# If reads_to_process not defined, set to 0 to process all reads
try: rtp = config['rtp'] 
except: rtp=0

try: genome_index = "global/"+config['idx']
except: genome_index = glob.glob("global/*.1.ht2")[0].replace(".1.ht2", "")

try:
  os.mkdir("fastp")
  os.mkdir("mapped")
  os.mkdir("counts")
except:
  None
########################## Rules ###########################

rule all:
    input:
        expand("counts/{sample}.counts", sample=SAMPLES),

if fastq_type == "single":
    rule run_fastp:
        input:
            "{sample}.fastq.gz",
        output:
            "fastp/{sample}.fastq.gz",
            "reports/{sample}.html"      
            #temp("mapped/{sample}.{params.fqs}"),
            #"reports/{sample}.html"
        params:
            #fqs = fastq_suffix,
            rtp = rtp,
        shell:
            """fastp\
            --reads_to_process {params.rtp}\
            --in1 {input[0]} --out1 {output[0]}\
            --html {output[1]}"""
    
    rule run_hisat2:
        input:
            "fastp/{sample}.fastq.gz",
        output:
            "mapped/{sample}.bam"
        params:
            idx = genome_index
        shell:
            'hisat2 -p 32 --max-intronlen 6000 -x {params.idx} -U {input} | sambamba view -S -f bam -o /dev/stdout /dev/stdin | sambamba sort --tmpdir="tmpmba" -t 32 -o {output} /dev/stdin'

    rule run_featureCounts:
        input:
            "mapped/{sample}.bam",
        output:
            "counts/{sample}.counts"
        params:
            idx = genome_index
        shell:
            "featureCounts -t exon,CDS -T 32 -a {params.idx}.gtf -o {output} {input}"
 

# Run fastp trimming on the fastq file
if fastq_type == "paired":
  rule run_fastp:
    input:
      expand("{sample}_1.fastq.gz", sample=SAMPLES, fqs=fastq_suffix),
      expand("{sample}_2.fastq.gz", sample=SAMPLES, fqs=fastq_suffix),
    output:
      "mapped/{sample}_1.fastp.gz",
      "mapped/{sample}_2.fastp.gz",
      "reports/{sample}.html"
    params:
      rtp = rtp,
      fqs = fastq_suffix
    shell:
      """mapped\
      --reads_to_process {params.rtp}\
      --in1 {input[0]}.{params.fqs} --out1 {output[0]}\
      --in2 {input[1]}.{params.fqs} --out2 {output[1]}\
      --html {output[2]}"""


if fastq_type == "paired":
  rule run_hisat2:
    input:
      "mapped/{sample}.fastp.gz",
    output:
      "mapped/{sample}.bam"
    params:
      hs2 = genome_index
    shell:
      """hisat2 -p 32 --max-intronlen 6000\
      -x {params.hs2} -U {input} | \
      sambamba view -S -f bam -o /dev/stdout /dev/stdin | \
      sambamba sort  --tmpdir="tmpmba" -t 32 -o {output} /dev/stdin
      """
