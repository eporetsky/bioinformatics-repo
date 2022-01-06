import gzip

"""
This script was written to extract only gene and exon sequence coordinates from GFF3 files.
Currently only extracts T001 transcripts but should be modified to take longest CDS or custom list.
 
###################### Converts from this format #######################
1	gramene	gene	44289	49837	.	+	.	ID=gene:Zm00001d027230;biotype=protein_coding;description=Zm00001d027230;gene_id=Zm00001d027230;logic_name=maker_gene
1	gramene	mRNA	44289	49837	.	+	.	ID=transcript:Zm00001d027230_T001;Parent=gene:Zm00001d027230;biotype=protein_coding;transcript_id=Zm00001d027230_T001
1	gramene	five_prime_UTR	44289	44350	.	+	.	Parent=transcript:Zm00001d027230_T001
1	gramene	exon	44289	44947	.	+	.	Parent=transcript:Zm00001d027230_T001;Name=Zm00001d027230_T001.exon1;constitutive=1;ensembl_end_phase=0;ensembl_phase=-1;exon_id=Zm00001d027230_T001.exon1;rank=1
1	gramene	CDS	44351	44947	.	+	0	ID=CDS:Zm00001d027230_P001;Parent=transcript:Zm00001d027230_T001;protein_id=Zm00001d027230_P001
1	gramene	exon	45666	45803	.	+	.	Parent=transcript:Zm00001d027230_T001;Name=Zm00001d027230_T001.exon2;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001d027230_T001.exon2;rank=2
1	gramene	CDS	45666	45803	.	+	0	ID=CDS:Zm00001d027230_P001;Parent=transcript:Zm00001d027230_T001;protein_id=Zm00001d027230_P001
1	gramene	exon	45888	46133	.	+	.	Parent=transcript:Zm00001d027230_T001;Name=Zm00001d027230_T001.exon3;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001d027230_T001.exon3;rank=3
1	gramene	CDS	45888	46133	.	+	0	ID=CDS:Zm00001d027230_P001;Parent=transcript:Zm00001d027230_T001;protein_id=Zm00001d027230_P001
1	gramene	exon	46229	46342	.	+	.	Parent=transcript:Zm00001d027230_T001;Name=Zm00001d027230_T001.exon4;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001d027230_T001.exon4;rank=4
1	gramene	CDS	46229	46342	.	+	0	ID=CDS:Zm00001d027230_P001;Parent=transcript:Zm00001d027230_T001;protein_id=Zm00001d027230_P001
1	gramene	exon	46451	46633	.	+	.	Parent=transcript:Zm00001d027230_T001;Name=Zm00001d027230_T001.exon5;constitutive=1;ensembl_end_phase=0;ensembl_phase=0;exon_id=Zm00001d027230_T001.exon5;rank=5
1	gramene	CDS	46451	46633	.	+	0	ID=CDS:Zm00001d027230_P001;Parent=transcript:Zm00001d027230_T001;protein_id=Zm00001d027230_P001
1	gramene	exon	47045	47262	.	+	.	Parent=transcript:Zm00001d027230_T001;Name=Zm00001d027230_T001.exon6;constitutive=1;ensembl_end_phase=2;ensembl_phase=0;exon_id=Zm00001d027230_T001.exon6;rank=6
1	gramene	CDS	47045	47262	.	+	0	ID=CDS:Zm00001d027230_P001;Parent=transcript:Zm00001d027230_T001;protein_id=Zm00001d027230_P001
1	gramene	CDS	47650	47995	.	+	1	ID=CDS:Zm00001d027230_P001;Parent=transcript:Zm00001d027230_T001;protein_id=Zm00001d027230_P001
1	gramene	exon	47650	48111	.	+	.	Parent=transcript:Zm00001d027230_T001;Name=Zm00001d027230_T001.exon7;constitutive=1;ensembl_end_phase=-1;ensembl_phase=2;exon_id=Zm00001d027230_T001.exon7;rank=7
1	gramene	three_prime_UTR	47996	48111	.	+	.	Parent=transcript:Zm00001d027230_T001
1	gramene	exon	48200	49247	.	+	.	Parent=transcript:Zm00001d027230_T001;Name=Zm00001d027230_T001.exon8;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001d027230_T001.exon8;rank=8
1	gramene	three_prime_UTR	48200	49247	.	+	.	Parent=transcript:Zm00001d027230_T001
1	gramene	exon	49330	49837	.	+	.	Parent=transcript:Zm00001d027230_T001;Name=Zm00001d027230_T001.exon9;constitutive=1;ensembl_end_phase=-1;ensembl_phase=-1;exon_id=Zm00001d027230_T001.exon9;rank=9
1	gramene	three_prime_UTR	49330	49837	.	+	.	Parent=transcript:Zm00001d027230_T001

############################ to this format ############################
1	gramene	gene	44289	49837	.	+	.	ID=Zm00001d027230
1	gramene	exon	44289	44947	.	+	.	Parent=Zm00001d027230_T001
1	gramene	exon	45666	45803	.	+	.	Parent=Zm00001d027230_T001
1	gramene	exon	45888	46133	.	+	.	Parent=Zm00001d027230_T001
1	gramene	exon	46229	46342	.	+	.	Parent=Zm00001d027230_T001
1	gramene	exon	46451	46633	.	+	.	Parent=Zm00001d027230_T001
1	gramene	exon	47045	47262	.	+	.	Parent=Zm00001d027230_T001
1	gramene	exon	47650	48111	.	+	.	Parent=Zm00001d027230_T001
1	gramene	exon	48200	49247	.	+	.	Parent=Zm00001d027230_T001
1	gramene	exon	49330	49837	.	+	.	Parent=Zm00001d027230_T001
"""


gff3_dict = {
'ZxPI566673': "Zx-PI566673-REFERENCE-YAN-1.0_Zx00001a.1.gff3.gz",
'B73': "Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3.gz",
'Mo17': "Zm-Mo17-REFERENCE-CAU-1.0_Zm00014a.1.gff3.gz",
'W22': "Zm-W22-REFERENCE-NRGENE-2.0_Zm00004b.1.gff3.gz",
'P39': "Zm-P39-REFERENCE-NAM-1.0_Zm00040a.1.gff3.gz",
'Oh43': "Zm-Oh43-REFERENCE-NAM-1.0_Zm00039a.1.gff3.gz",
'NC358': "Zm-NC358-REFERENCE-NAM-1.0_Zm00037a.1.gff3.gz",
'NC350': "Zm-NC350-REFERENCE-NAM-1.0_Zm00036a.1.gff3.gz",
'Mo18W': "Zm-Mo18W-REFERENCE-NAM-1.0_Zm00034a.1.gff3.gz",
'M37W': "Zm-M37W-REFERENCE-NAM-1.0_Zm00032a.1.gff3.gz",
'Ky21': "Zm-Ky21-REFERENCE-NAM-1.0_Zm00031a.1.gff3.gz",
'Ki3': "Zm-Ki3-REFERENCE-NAM-1.0_Zm00029a.1.gff3.gz",
'Ki11': "Zm-Ki11-REFERENCE-NAM-1.0_Zm00030a.1.gff3.gz",
'Il14H': "Zm-Il14H-REFERENCE-NAM-1.0_Zm00028a.1.gff3.gz",
'HP301': "Zm-HP301-REFERENCE-NAM-1.0_Zm00027a.1.gff3.gz",
'CML69': "Zm-CML69-REFERENCE-NAM-1.0_Zm00020a.1.gff3.gz",
'CML52': "Zm-CML52-REFERENCE-NAM-1.0_Zm00019a.1.gff3.gz",
'CML333': "Zm-CML333-REFERENCE-NAM-1.0_Zm00026a.1.gff3.gz",
'CML277': "Zm-CML277-REFERENCE-NAM-1.0_Zm00024a.1.gff3.gz",
'CML247': "Zm-CML247-REFERENCE-NAM-1.0_Zm00023a.1.gff3.gz",
'CML228': "Zm-CML228-REFERENCE-NAM-1.0_Zm00022a.1.gff3.gz",
'B97': "Zm-B97-REFERENCE-NAM-1.0_Zm00018a.1.gff3.gz",
'Tzi8': "Zm-Tzi8-REFERENCE-NAM-1.0_Zm00042a.1.gff3.gz",
'Tx303': "Zm-Tx303-REFERENCE-NAM-1.0_Zm00041a.1.gff3.gz",
'Oh7B': "Zm-Oh7B-REFERENCE-NAM-1.0_Zm00038a.1.gff3.gz",
'Ms71': "Zm-Ms71-REFERENCE-NAM-1.0_Zm00035a.1.gff3.gz",
'M162W': "Zm-M162W-REFERENCE-NAM-1.0_Zm00033a.1.gff3.gz",
'CML322': "Zm-CML322-REFERENCE-NAM-1.0_Zm00025a.1.gff3.gz",
'CML103': "Zm-CML103-REFERENCE-NAM-1.0_Zm00021a.1.gff3.gz"
}

write_mRNA = False
mRNA = "001"
for key, value in gff3_dict.items():
    with gzip.open("../GFF3/"+value, "rt") as handle:
        #with gzip.open("gff3/"+key+".gff3.gz","wt",compresslevel=9) as fixed:
        with open("gff3/"+key+".gff3","w+") as fixed:
            for row in handle:
                if "\t" not in row or row[0] == "#":
                    continue

                row = row.split("\t")
                desc = row[-1].split(";")

                if row[2] == "gene":
                    gene = desc[0].split("=")[1]
                    if gene.startswith("gene:"):
                        gene = gene[5:]
                    row[-1] = "=".join(["ID",gene])
                    fixed.write("\t".join(row)+"\n")
                    continue

                if row[2] == "mRNA":
                    mRNA = desc[0].split("=")[1]
                    if mRNA.startswith("transcript:"):
                        mRNA = mRNA[11:]
                    row[-1] = "=".join(["ID",mRNA])
                    if write_mRNA:
                        fixed.write("\t".join(row)+"\n")
                    continue

                if mRNA[-3:] != "001":
                    continue

                if row[2] == "exon":
                    row[-1] = "=".join(["Parent",mRNA])
                    fixed.write("\t".join(row)+"\n")
                    continue
