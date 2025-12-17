import pandas as pd

df = pd.read_table("0.human_gene_info_from_Ensembl.txt", sep = "\t")

tmp = list(set(list(df['Gene name'])).intersection(list(df['Gene Synonym'])))
tmp = [i for i in tmp if str(i) != 'nan'] # remove NaN

filter_ = []
for i in tmp:
    tmp1 = list(set(df[df['Gene name'] == i]['Gene stable ID']))
    tmp2 = list(set(df[df['Gene Synonym'] == i]['Gene stable ID']))
    check = [x for x in tmp2 if x not in tmp1]
    if check != []:
        filter_.append(i)

df.loc[df['Gene Synonym'].isin(filter_), 'Gene Synonym'] = ''


df.to_csv('0.human_gene_info_from_Ensembl.mod', sep='\t', index=False)

