# LRRfinder
Following the installation process here: https://github.com/cgottin/LRRprofiler\
```
gzip --decompress --keep  protein_seqs/Zmays_493_RefGen_V4.protein_primaryTranscriptOnly.locus.fa.gz
singularity run LRRprofiler.sif --in_proteome protein_seqs/Zmays_493_RefGen_V4.protein_primaryTranscriptOnly.locus.fa --name B73v4 --nobuild
```

```
To install phobius, tmhmm and signalp:

conda install -c predector phobius
Register: https://software.sbc.su.se/cgi-bin/request.cgi?project=phobius
phobius-register phobius101_linux.tar.gz

conda install -c predector signalp5
Register: https://services.healthtech.dtu.dk/service.php?SignalP-5.0
signalp5-register signalp-5.0b.Linux.tar.gz


conda install -c predector tmhmm
# Register: https://services.healthtech.dtu.dk/service.php?TMHMM-2.0
tmhmm2-register tmhmm-2.0c.Linux.tar.gz
```

# InterProScan

### 1. Register SignalP-4.1
```
# You will first need to register and download SignalP-4.1 (no need to conda install)
# https://services.healthtech.dtu.dk/service.php?SignalP-4.1

# Download interperoscan standalone version
# This version requires phobius, tmhmm and signalp. Since you had to download the standalone versions when registring,
# it seems easier to just copy the 3 packages to the existing data/ folder and create a packages/ folder in it.
```

### 2. Then edit interproscan.properties (apologies for not uploading the proporties file itself, not sure what the )

```

#signalp
# Note: SignalP binary not distributed with InterProScan 5, please install separately e.g. in ${bin.directory}/signalp/4.1/signalp
binary.signalp.path=${bin.directory}/signalp/4.1/signalp
signalp.perl.library.dir=${bin.directory}/signalp/4.1/lib

#TMHMM 2.0
# Note: TMHMM binary not distributed with InterProScan 5, please install separately e.g. in ${bin.directory}/tmhmm/2.0c/decodeanhmm
binary.tmhmm.path=${bin.directory}/tmhmm/2.0c/decodeanhmm

#PHOBIUS
# Note: Phobius binary not distributed with InterProScan 5, please install separately e.g. in ${bin.directory}/phobius/1.01/phobius.pl
binary.phobius.pl.path=${bin.directory}/phobius/1.01/phobius.pl

#TMHMM 2.0
# Note: TMHMM model files not distributed with InterProScan 5, please install separately e.g. in data/tmhmm/2.0/TMHMM2.0.model
tmhmm.model.path=${data.directory}/tmhmm/2.0c/TMHMM2.0c.model
```

### 3. Clean fasta file for InterProScan

```
# decompress all fasta files in a folder. Assumes ".fa" suffix, edit file otherwise.
# This script removes id descriptions and any * from the seqs. Check output before proceeding. 
python3 clean_fasta.py all
```

### 4. Run InterProScan

```
# Test InterProScan:
interproscan.sh -i test_all_appl.fasta -f tsv
# 100% done:  InterProScan analyses completed

To Run InterProScan:
./interproscan.sh -i proteomes/protein_seqs.clean.fa -f tsv,GFF3 -cpu 32 -dp -dra
```