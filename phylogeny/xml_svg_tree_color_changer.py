
name_list=[]
annotated=0

f = open('Aligned/GRMZM2G451007_T001_tree','r')
w = open('Aligned/GRMZM2G451007_T001_colored_tree.svg','w')

count = 0 
for line in f:
    if 313 > count > 5:
        name = line[line.find(">")+1:line.find("</text>")]
        if name in responders:
            w.write(line[:line.find("rgb(")] + "rgb(255,0,0)" + line[line.find("rgb(")+10:])
            continue
    w.write(line)            
    count+=1
w.close()