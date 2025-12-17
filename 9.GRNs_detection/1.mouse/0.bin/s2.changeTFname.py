
dic = []
for line in open('allTFs_mm.txt'):
    dic.append(line.strip())

flag = 0
change= {}
for line  in open('0.mouse_gene_info_from_Ensembl.mod'):
    flag += 1
    if flag >= 2:
        line = line.strip().split('\t')
        if len(line) == 1:
            line = [line[0], line[0], line[0]]
        if len(line) == 2:
            line.append(line[1])
        if len(line) ==3:
            line = [i if i != "" else line[0] for i in line]

        try:
            change[line[1]].append(line[0])
        except:
            change[line[1]] = [line[0]]

        try:
            change[line[2]].append(line[0])
        except:
            change[line[2]] = [line[0]]
out = []
write_ = open("TFs_Ensembl.txt", 'w')
for i in dic:
    try:
        out.extend(list(set(change[i])))
    except:
        out.append(i)
write_.write('\n'.join(list(set(out))))
write_.close()
