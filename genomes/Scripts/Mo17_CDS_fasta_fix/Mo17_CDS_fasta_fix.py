from Bio import SeqIO
from Bio.SeqIO import SeqRecord

"""
Link to Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.cds.fa.gz:
https://ftp.maizegdb.org/MaizeGDB/FTP/Zm-Mo17-REFERENCE-CAU-1.0/

Below is an example how the CDS fasta file is annotated:
>Zm00014a008050_T001 gene=Zm00014a008050 CDS=57-1985
The python script below retrives the annotated CDS sequence

print(count_ATG) returns [45482, 1048], meaning:
1) 45k genes have ATG at annotated CDS coordinates
2) 1k don't have an ATG there
print(count_stop) returns [46527, 3], meaning:
All but 3 genes don't have a stop codon at the annotated site

All annotated CDS sequences were were "fixed", regardless of ATG or stop presence
Hope this might be useful to other people
"""

count_ATG=[0,0]
count_stop=[0,0]

gene_dict = dict()
records = list(SeqIO.parse("Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.CDS.fa", "fasta"))
for record in records:
	description=record.description
	name = record.name
	coords = description.split(" ")[-1][4:].split("-")
	seq = record.seq # protein sequence
	cds = seq[int(coords[0])-1:int(coords[1])]
	gene_dict[description] = SeqRecord(seq=cds, id=name, description=description)

	if cds[0:3]=="ATG":
		count_ATG[0]+=1
	else:
		count_ATG[1]+=1
	if cds[-3:] in ["TAA","TAG","TGA"]:
		count_stop[0]+=1
	else:
		count_stop[1]+=1
print("Count of start codons:", count_ATG)
print("Count of stop codons:", count_stop)

with open("Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.CDS.fixed.fa", "w") as handle:
	SeqIO.write(gene_dict.values(), handle, "fasta")
