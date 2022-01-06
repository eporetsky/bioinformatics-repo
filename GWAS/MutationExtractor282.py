# https://pyvcf.readthedocs.io/en/latest/INTRO.html
 
from collections import Counter
from gff3 import Gff3
import vcf
import csv

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

gene_list = ["GRMZM2G168629","GRMZM2G471304","GRMZM2G168888","GRMZM2G168984","GRMZM2G169009","GRMZM2G169020","GRMZM2G028980","GRMZM2G433601","GRMZM2G132854","GRMZM2G433609","GRMZM2G132862","GRMZM2G132866","GRMZM2G401831","GRMZM2G401869","GRMZM2G101818","GRMZM2G101872","GRMZM2G101916","GRMZM2G101920","GRMZM2G434191","GRMZM2G118637","GRMZM2G419844","GRMZM2G419891","GRMZM2G119640","GRMZM2G419994","GRMZM5G834952","GRMZM2G420004","GRMZM5G861442","GRMZM2G119932","GRMZM2G420010","GRMZM2G111293","GRMZM2G111304","GRMZM2G111319","GRMZM2G410963","GRMZM2G012262","GRMZM2G012143","GRMZM2G009871","GRMZM2G052740","GRMZM2G052720","GRMZM2G052696","GRMZM2G351027","GRMZM2G003157","GRMZM2G303374","GRMZM5G888445","GRMZM2G450937","GRMZM2G450937","GRMZM2G451007","GRMZM2G451007","GRMZM2G153745","GRMZM2G153797","GRMZM5G851807","GRMZM2G153899","GRMZM2G058004","GRMZM2G166162","GRMZM2G467703","GRMZM2G467682","GRMZM2G166147","GRMZM2G166145","GRMZM2G166044","GRMZM2G166041","GRMZM2G166035","GRMZM2G166032","GRMZM2G165917","GRMZM2G155370","GRMZM2G155384","GRMZM2G117346","GRMZM2G117441","GRMZM2G117459","GRMZM2G165931","GRMZM2G165966","GRMZM2G065494","GRMZM2G065451","GRMZM2G092627","GRMZM2G139341","GRMZM2G139372","GRMZM2G123511","GRMZM5G846811","GRMZM2G318515","GRMZM2G389768","GRMZM2G090156","GRMZM2G042602","GRMZM2G042664","GRMZM2G061078","GRMZM2G061023","GRMZM2G061016","GRMZM2G060669","GRMZM2G060567","GRMZM2G060451","GRMZM2G362021","GRMZM2G423137","GRMZM2G423129","GRMZM2G121269","GRMZM2G174669","GRMZM2G174667","GRMZM2G174596","GRMZM2G174568","GRMZM2G174481","GRMZM2G133937","GRMZM2G133958","GRMZM2G133958","GRMZM2G133988","GRMZM2G435219","GRMZM5G867322","GRMZM2G435255","GRMZM5G871727","GRMZM5G882631","GRMZM5G851744","GRMZM5G866522","GRMZM5G848897","GRMZM2G178100","GRMZM5G879778","GRMZM2G382171","GRMZM5G818155","GRMZM2G027726","GRMZM2G027640","GRMZM2G027209","GRMZM2G027187","GRMZM2G027043","GRMZM2G319875","GRMZM2G419209","GRMZM2G010044"]
#gene_list = ["GRMZM2G451007", "GRMZM2G303374"]

def coordinate_list_generator(gene_list):
    isoforms = ["_T001","_T002","_T003","_T004","_T005","_T006","_T007","_T008","_T009","_T010","_T011","_T012","_T013","_T014","_T015","_T016","_T017","_T018","_T019","_T020","_T021","_T022","_T023","_T024","_T025","_T026","_T027","_T028","_T029","_T030"]
    with open('sequence.gff3','r') as tsv:
        AoA = [line.strip().split('\t') for line in tsv]

    coordinate_dict = {}
    gene_direction = {}
    count=0
    for row in range(2,len(AoA)-1):

        i = AoA[row]
        if i[2]=="CDS":
            gene_loc = i[8].find("gene=")
            if i[8][gene_loc+5:gene_loc+8]=="GRM":
                if i[8][gene_loc+5:gene_loc+18] in gene_list:
                    gene_name = i[8][gene_loc+5:gene_loc+18]+isoforms[count]
                    coordinate_dict.setdefault(gene_name, []).append([i[3], i[4]])
                    gene_direction[gene_name]=i[6]
                    if AoA[row+1][2]=="gene":
                        count=0
                        continue
                    elif AoA[row+1][2]=="mRNA":
                        count+=1
            else:
                if i[8][gene_loc+5:gene_loc+21] in gene_list:
                    gene_name = i[8][gene_loc+5:gene_loc+21]+isoforms[count]
                    coordinate_dict.setdefault(gene_name, []).append([i[3], i[4]])
                    gene_direction[gene_name]=i[6]
                    if AoA[row+1][2]=="gene":
                        count=0
                        continue
                    elif AoA[row+1][2]=="mRNA":
                        count+=1
    return (coordinate_dict, gene_direction)
    

def create_gene_fasta(coordinate_dict, gene_direction, chromosome_seq):
    ### Saves the sequences for every gene in the provided list as a fasta file
    for gene_name, coord_list in coordinate_dict.iteritems():
        cds_seq_dict = {}
        for line in generate_line_list():
            cds_seq_dict[line] = ""
        for exon_coord in coord_list:
            exon_seq = create_exon_fasta(exon_coord, chromosome_seq, gene_direction[gene_name])
            cds_seq_dict = dict(Counter(cds_seq_dict) + Counter(exon_seq))
        ### Below takes a dictionary of sequences for a gene and saves it as a fasta file
        fasta_filename = gene_name+".fasta"
        output_handle = open(fasta_filename, "w")
        for line_name, seq in cds_seq_dict.iteritems():
            output_handle.write(">%s \n%s\n" % (line_name, seq))
        output_handle.close()

#create_gene_fasta(coordinate_dict, chromosome_seq)        

def create_exon_fasta(exon_coord, chromosome_seq, direction):
    ### This is the main function that takes one exon-cds, finds all the SNPs within it
    ### and creates a single dictionary object with an entry for every corrected line
    ##### Need to find a way to correct for indels!!!
    
    left, right = int(exon_coord[0])-1, int(exon_coord[1])
    exon_seq = chromosome_seq[left:right].seq.__str__()

    
    seq_dict = {} # Stores the sequences for every variety/line in the vcf file
    for line in generate_line_list():
        seq_dict[line] = exon_seq

    vcf_reader = vcf.Reader(open('282_FAC.vcf', 'r'))
    for snp in vcf_reader:
        if left <= snp.POS <= right:
            snp_left, snp_right = snp.POS-1, snp.POS+len(snp.REF)-1
            ref, alt = snp.REF, snp.ALT[0]
            for line in snp:
                if line['GT'] == "1/1":
                    temp_seq = seq_dict[line.sample]
                    seq_dict[line.sample] = temp_seq[:snp_left-left]+alt.__str__()+temp_seq[snp_right-left:]
                else:
                    continue
        else:
            continue
    if direction=="-":
        rev_seq_dict = {k: reverse_complement(v) for k, v in seq_dict.items()}
        return rev_seq_dict
    return seq_dict
    
    

def generate_line_list(vcf_file='282_receptors.vcf'):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    vcf_record = next(vcf_reader)
    line_list = []
    for line in vcf_record:
        line_list.append(line.sample)
    return line_list

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])
    
    
from Bio import SeqIO
seq_list = []
name_list = ["names"]
with open("Done/282_FAC_Gene_Seqs_T001/GRMZM2G451007_T001.fasta", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        seq_list.append(str(record.seq.translate()))
        name_list.append(str(record.name))

seq_len = len(seq_list[0])
mut_list = [name_list]
count = 0
for i in range(seq_len):
    if len(set(list(zip(*seq_list))[i])) > 1:
        mut_list.append([str(i)])
        for aa in zip(*seq_list)[i]:
            mut_list[count+1].append(aa)
        count += 1
#for aa in str(record.seq.translate()):
#print aa
#break
import csv
mut_list_t = zip(*mut_list)
with open("output.csv", "wb") as f:
    writer = csv.writer(f)
    writer.writerows(mut_list_t)