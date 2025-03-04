import loompy
import scanpy as sc
import pandas
import numpy
import scipy
from scipy import io
import random
import numpy as np

adata_ = sc.read_loom('l5_all.loom')

#### only include protein coding genes ####
pep = []
for line in open('GeneID_name.txt'):
    line =line.strip().split('\t')
    if line[-1] == 'protein_coding':
        pep.append(line[1])

genes_to_keep = adata_.var_names.isin(pep)
adata = adata_[:, genes_to_keep].copy()


# Move embeddings info to the right place and right format
x = pandas.Series.to_numpy(adata.obs['_X'])
y = pandas.Series.to_numpy(adata.obs['_Y'])
xy = numpy.stack((x,y)).transpose().reshape(-1,2)
adata.obsm['X_test'] = xy

pc1 = pandas.Series.to_numpy(adata.obs['_PC1'])
pc2 = pandas.Series.to_numpy(adata.obs['_PC2'])
pc = numpy.stack((pc1,pc2)).transpose().reshape(-1,2)
adata.obsm['PCA'] = pc

t1 = pandas.Series.to_numpy(adata.obs['_tSNE1'])
t2 = pandas.Series.to_numpy(adata.obs['_tSNE2'])
tsne= numpy.stack((t1,t2)).transpose().reshape(-1,2)
adata.obsm['tSNE'] = tsne

# Only include necessary metadata:
adata.obs = adata.obs[[ 'Age', 'AnalysisPool', 'AnalysisProject', 'Class', 'ClusterName', 'Clusters', 'Description', 'DonorID', 'MitoRiboRatio', 'Neurotransmitter', 'TaxonomyRank1', 'TaxonomyRank2', 'TaxonomyRank3', 'TaxonomyRank4', 'TaxonomySymbol', 'Tissue' ]] 

# Change the matrix format
adata.X = scipy.sparse.csc_matrix(adata.X)

# Make variable and observation names unique
adata.var_names_make_unique()
adata.obs_names_make_unique()


# subset only brain, Pons also excluded
subset = adata[adata.obs.Tissue.isin(["Amygd", "CA1", "CB", "Ctx1", "Ctx1.5", "Ctx2", "Ctx3", "DentGyr", "HC", "Hypoth", "MBd", "MBv", "Medulla", "OB", "SScortex", "StriatDor", "StriatVent", "Thal"])]

# gene need to be deleted
deleteID = []
for line in open('IDdelete'):
    if not line.startswith('#'):
        deleteID.append(line.strip())
subset = subset[:, ~subset.var_names.isin(deleteID)]
# Write h5ad file
subset.write('Mmus.h5ad')


# Wirte data for Seurat processing
with open('matrix_files/barcodes.tsv', 'w') as f:
    for item in subset.obs_names:
        f.write(item + '\n')

with open('matrix_files/features.tsv', 'w') as f:
    for item in ['\t'.join([x, x, "Gene Expression"]) for x in subset.var_names]:
        f.write(item + '\n')

io.mmwrite('matrix_files/matrix', subset.X.T)
subset.obs.to_csv('metadata100k.csv')

