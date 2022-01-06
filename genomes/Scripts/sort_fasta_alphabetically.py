import glob

from Bio import SeqIO


"""
Finds all .fa files in a folder and sorts the ids alphabetically
"""

fasta_list = glob.glob('*.fa') # returns a list of all .fa files in the folder
for fasta in fasta_list: 
	print(str("Sorting "+fasta))
	# convert fasta to dictionary
	record_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))	
	# return list of keys and sort alphabetically
	sorted_list = sorted(record_dict.keys())
    # write sequences to a file in the order of the sorted_list
	with open(str(fasta[:-3]+"_sorted.fa"), "w") as f:
		for gene in sorted_list:
			f.write(str(">"+gene+"\n"))
			f.write(str(record_dict[gene].seq+"\n"))
