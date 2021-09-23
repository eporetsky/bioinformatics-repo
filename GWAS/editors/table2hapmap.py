# This is how the original file looks:
# Pos	B73v4_Ref	line1	line2	line3	
# 1:2400	A	N	N	N	

# That's how it should look:
# rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode	line1	line2	line3	
# ss1_3498	G/A	1	3498	+	NA	NA	NA	NA	NA	NA	N	N	N

import csv

input_file = open("input_hapmap.txt", "r")
with open('output_hapmap.txt', 'w', newline='') as f_output:
    
    tsv_output = csv.writer(f_output, delimiter='\t')
    n=0

    for aline in input_file:
        n+=1
        if n<2:
            spl = aline.split("\t")
            new_line = ["rs#","allele","chrom","pos","strand","assembly#","center","protLSID","assayLSID","panelLSID","QCcode"]
            new_line.extend(spl[2:])
            new_line[-1] = new_line[-1][:-1]
            tsv_output.writerow(new_line)
    
        if n>1:
            spl = aline.split("\t")
            alleles = set(spl[1:-1])
            try:
                alleles.remove("N")
            except:
                continue
            alleles.remove(spl[1])
            try:
                alleles = str(spl[1]+"/"+(list(alleles)[0]))
            except:
                alleles = str(spl[1]+"/N")
            chrom_pos = spl[0].split(":")
    
            new_line = [spl[0],alleles, chrom_pos[0], chrom_pos[1], "+","NA","NA","NA","NA","NA","NA"]
            new_line.extend(spl[2:])
            new_line[-1] = new_line[-1][:-1]
            
            tsv_output.writerow(new_line

input_file.close()


