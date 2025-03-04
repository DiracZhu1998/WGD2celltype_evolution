from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            CellTypeTriangles, ParalogSubstitutions, GeneTriangles)
import samalg
import pandas as pd
import scanpy as sc
import numpy as np
import pickle
from igraph import Graph
import sys, re, os


atlas = sys.argv[1]
label = sys.argv[2]
outdir = sys.argv[3]
os.makedirs(outdir, exist_ok=True)

with open(atlas, 'rb') as f:
    sm = pickle.load(f)

# add refined annotation to SAMap object for comparison

info1 = pd.read_csv('Refined_annotation/Pmar.iter_clustering.ann.count_mean.txt',sep = '\t',names = ['cell','iter.clustering'], index_col = 0)
info2 = pd.read_csv('Refined_annotation/Pvit.iter_clustering.ann.count_mean.txt',sep = '\t',names = ['cell','iter.clustering'], index_col = 0)
info3 = pd.read_csv('Refined_annotation/Mmus.iter_clustering.ann.count_mean.txt',sep = '\t',names = ['cell','iter.clustering'], index_col = 0)
info4 = pd.read_csv('Refined_annotation/Hsap.iter_clustering.ann.count_mean.txt',sep = '\t',names = ['cell','iter.clustering'], index_col = 0)

sm.sams['pm'].adata.obs['iter.clustering'] = info1['iter.clustering']
sm.sams['pv'].adata.obs['iter.clustering'] = info2['iter.clustering']
sm.sams['mm'].adata.obs['iter.clustering'] = info3['iter.clustering']
sm.sams['hs'].adata.obs['iter.clustering'] = info4['iter.clustering']


Pmar = pd.read_csv("Refined_annotation/Pmar_family_annotation.txt", sep = "\t")
Mmus = pd.read_csv("Refined_annotation/Mmus_family_annotation.txt", sep = "\t")
Pvit = pd.read_csv("Refined_annotation/Pvit_family_annotation.txt", sep = "\t")
Hsap = pd.read_csv("Refined_annotation/Hsap_family_annotation.txt", sep = "\t")

result = pd.merge(sm.sams['pm'].adata.obs, Pmar, left_on='iter.clustering', right_on='Cluster', how='left')
result.index = sm.sams['pm'].adata.obs.index
sm.sams['pm'].adata.obs = result

result = pd.merge(sm.sams['mm'].adata.obs, Mmus, left_on='iter.clustering', right_on='Cluster', how='left')
result.index = sm.sams['mm'].adata.obs.index
sm.sams['mm'].adata.obs = result

result = pd.merge(sm.sams['pv'].adata.obs, Pvit, left_on='iter.clustering', right_on='Cluster', how='left')
result.index = sm.sams['pv'].adata.obs.index
sm.sams['pv'].adata.obs = result

result = pd.merge(sm.sams['hs'].adata.obs, Hsap, left_on='iter.clustering', right_on='Cluster', how='left')
result.index = sm.sams['hs'].adata.obs.index
sm.sams['hs'].adata.obs = result


keys = {key: 'Refined subtype' for key in sm.sams.keys()}
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)

#### generate linkage data for R chord plot
A = MappingTable
xx=A.values.copy()
x,y = xx.nonzero()
z=xx[x,y]
x,y = A.index[x],A.columns[y]
links=pd.DataFrame(data=np.array([x,y,z]).T,columns=['source','target','value'])
links = links[links['value'] >= 0.1] # mapping score >= 0.1

""" previous test, group cells based on linkage
# use igraph to link cell types with linkage into same groups
#g = Graph.TupleList(links[['source', 'target']].itertuples(index=False), directed=False, weights=False)
#components = g.connected_components()
#membership = components.membership

#groups = {}
#for idx, vertex in enumerate(g.vs):
#    comp_id = membership[idx]
#    groups[vertex['name']] = str(comp_id)
"""
#links['groups'] = links['source'].map(groups).astype('str')
# filter clusters with less than 100 cells and filter Mixture clusters
# neural and non-neural atlas, few neurons for example could be misclustered into glial clusters but identified in my iterative clustering

filter_list = []
def filtering(tmp, label):
    ID1 = (tmp[tmp <= 100].index.tolist())
    ID1 = [label + '_' + i for i in ID1]
    ID2 = [label + '_' + i for i in tmp.index.tolist() if 'Mixture' in i]
    IDs = list(set(ID1 + ID2))
    return IDs

filter_list.extend(filtering(sm.sams['pm'].adata.obs['Refined subtype'].value_counts(), 'pm'))
filter_list.extend(filtering(sm.sams['pv'].adata.obs['Refined subtype'].value_counts(), 'pv'))
filter_list.extend(filtering(sm.sams['mm'].adata.obs['Refined subtype'].value_counts(), 'mm'))
filter_list.extend(filtering(sm.sams['hs'].adata.obs['Refined subtype'].value_counts(), 'hs'))

filtered_links = links[~links['source'].isin(filter_list)]
filtered_links = filtered_links[~filtered_links['target'].isin(filter_list)]
links.to_csv(f'{outdir}/{label}.MappingTables.links.csv', index=False)
filtered_links.to_csv(f'{outdir}/{label}.MappingTables.links_filtered.csv', index=False)

sc.tl.leiden(sm.samap.adata, resolution=1.0, key_added = "SAMap.integrated_leiden.1.0")
sc.tl.leiden(sm.samap.adata, resolution=1.5, key_added = "SAMap.integrated_leiden.1.5")
sc.tl.leiden(sm.samap.adata, resolution=2.0, key_added = "SAMap.integrated_leiden.2.0")

umap = pd.DataFrame(sm.samap.adata.obsm['X_umap'].copy(), columns=['UMAP1', 'UMAP2'], index=sm.samap.adata.obs.index)

umap.to_csv(f'{outdir}/{label}.UMAP.csv', index = True)
sm.samap.adata.obs.to_csv(f'{outdir}/{label}.meta.csv', index=True)
MappingTable.to_csv( f'{outdir}/{label}.MappingTable.csv' , index=True)

# save dataset
import pickle
with open(f'{outdir}/{label}.pkl', 'wb') as f:  # open a text file
    pickle.dump(sm, f)
f.close()

# generate gene pairs supporting cell type pairs:
#gpf = GenePairFinder(sm,keys=keys)
#gene_pairs = gpf.find_all(align_thr=0.1)
#gene_pairs.to_csv(f'{outdir}/{label}.samap.gene_pairs.csv', index = False)
