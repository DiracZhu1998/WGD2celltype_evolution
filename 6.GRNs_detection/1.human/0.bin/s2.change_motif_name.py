

flag = 0
change= {}
for line  in open('0.human_gene_info_from_Ensembl.mod'):
    flag += 1
    if flag >= 2:
        line = line.strip().split('\t')
        if len(line) == 1:
            line = [line[0], line[0], line[0]]
        if len(line) == 2:
            line.append(line[1])
        if len(line) ==3:
            line = [i if i != "" else line[0] for i in line]

        #print(line)


        try:
            change[line[1]].append(line[0])
        except:
            change[line[1]] = [line[0]]

        try:
            change[line[2]].append(line[0])
        except:
            change[line[2]] = [line[0]]


for line in open('motifs-v9-nr.hgnc-m0.001-o0.0.tbl'):
    if line.startswith('#'):
        print(line.strip())
    else:
        line = line.strip().split('\t')
        try:
            gene = list(set(change[line[5]]))
        except:
            gene = [line[5]]

        for i in gene:
            line[5] = i
            print('\t'.join(line))


