from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import glob
fasta_list = glob.glob('./*.fa') # returns a list of all .fa files in the folder

for fasta in fasta_list: 
	# genome_name is a maize specific containing genome number and name (i.e. 14_Mo17)
	genome_name = str(fasta.split("-")[4][9:11]+"_"+fasta.split("-")[1]+".fa") 
	records = list(SeqIO.parse(fasta, "fasta")) # extract all records to list
	gene_dict = dict()
	# iterate over records and check for matching gene ids
	for record in records:
		name = (record.name)[0:14] # gene id excluding the isoform
		seq = record.seq # protein sequence
		ln = len(seq) # length of protein sequence
		if name in gene_dict:
			# if the new isoform is longer than previous one, replace it in dictionary
			if ln > len(seq):
				gene_dict[name] = SeqRecord(seq=seq, id=name, description="")
		else:
			gene_dict[name] = SeqRecord(seq=seq, id=name, description="")

	# write dictionary of SeqRecors to fasta file
	with open(genome_name, "w") as handle:
		SeqIO.write(gene_dict.values(), handle, "fasta") 