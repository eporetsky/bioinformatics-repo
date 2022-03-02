import sys
import glob
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# This script can be used to clean all fasta sequences in a folder.
# Alternatively, use primary2locus.ipynb to generate a clean fasta

# Use args, only a single one for now, to clean a single selected fasta file
if len(sys.argv) == 2:
	file_name = sys.argv[1]
else:
	file_name = "all"

if file_name == "all":
	fasta_list = glob.glob('*.fa') # returns a list of all .fa files in the folder
else:
	fasta_list = [file_name]

for fasta in fasta_list: 
	print(str("Processing "+fasta))
	records = list(SeqIO.parse(fasta, "fasta")) # extract all records to list
	gene_dict = dict()
	# iterate over records and check for matching gene ids
	for record in records:
		name = record.name.split("_")[0] # gene id excluding the isoform. Might not work for all genomes
		seq = record.seq # protein sequence
		# InterProScan can't handle * in the fasta file. Might as well remove
		# Word of caution - some times fastsa sequences have * in the middle of the AA seq, not just the end fo seq
		# I don't know if this is an error, but the following script will remove all * from the fasta file, indiscriminately
		gene_dict[name] = SeqRecord(seq=seq.replace("*",""), id=name, description="")

	# write dictionary of SeqRecors to fasta file
	with open(fasta.replace(".fa","")+".clean.fa", "w") as handle:
		SeqIO.write(gene_dict.values(), handle, "fasta") 
