import anndata, random
import scanpy as sc
from scipy import io
import pandas as pd
import numpy as np

adata = sc.read_h5ad("human_adult_GRCh38-3.0.0.h5ad")


#### only include protein coding genes ####
pep = []
for line in open('GeneID_name.txt'):
    line =line.strip().split('\t')
    if line[-1] == 'protein_coding':
        pep.append(line[0])

# only brain
subset = adata[~adata.obs.roi.isin(["SpC"])] 
genes_to_keep = subset.var_names.isin(pep)
subset = subset[:, genes_to_keep].copy()
sc.pp.filter_cells(subset, min_counts = 400)


# read Supp. Table and more information
super_info = pd.read_csv("science.add7046_table_s3.csv")
subset.obs['supercluster'] = pd.NA
subset.obs['cellname'] = pd.NA

for x in super_info['Cluster ID']:
    # Extract cluster_name and supercluster_name for this Cluster.ID
    cluster_name = super_info.loc[super_info['Cluster ID'] == x, 'Cluster name'].values[0]
    supercluster_name = super_info.loc[super_info['Cluster ID'] == x, 'Supercluster'].values[0]
    # Update corresponding cells in adata.obs
    subset.obs.loc[subset.obs['cluster_id'] == x, 'cellname'] = cluster_name
    subset.obs.loc[subset.obs['cluster_id'] == x, 'supercluster'] = supercluster_name


subset.obs['Region'] = [i.split(' - ')[0] for i in subset.obs['dissection'].tolist() ]
subset.obs['DonorID'] = subset.obs['donor_id']
# subset 100k neurons + 100k glail cells
list_ = ['Ependymal', 'Oligodendrocyte precursor', 'Microglia', 'Fibroblast',
        'Committed oligodendrocyte precursor', 'Bergmann glia', 'Astrocyte',
        'Oligodendrocyte', 'Vascular', 'Choroid plexus']
neurons_ID = subset.obs[~subset.obs['supercluster'].isin(list_)].sample(100000, random_state=42).index.tolist()
non_neurons_ID = subset.obs[subset.obs['supercluster'].isin(list_)].sample(100000, random_state=42).index.tolist()

subset = subset[subset.obs.index.isin(list(neurons_ID+non_neurons_ID))].copy()
subset.X = subset.X.astype(np.float64) # convert from int to float64
subset.write('Hsap.subset.100k_neurons_100k_others.v2.h5ad', compression = 'gzip')


with open('matrix_files/barcodes.tsv', 'w') as f:
    for item in subset.obs_names:
        f.write(item + '\n')

with open('matrix_files/features.tsv', 'w') as f:
    for item in ['\t'.join([x, x, "Gene Expression"]) for x in subset.var_names]:
        f.write(item + '\n')

io.mmwrite('matrix_files/matrix', subset.X.T)

subset.obs.to_csv('metadata.csv')
