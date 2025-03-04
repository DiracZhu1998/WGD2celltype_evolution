import sys
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import nsforest as ns
from nsforest import utils
from nsforest import preprocessing as pp
from nsforest import nsforesting
from nsforest import evaluating as ev
from nsforest import plotting as pl

adata_ = sys.argv[1]
cluster_header = sys.argv[2] # cluster name in annotation
output_folder = sys.argv[3]
label = sys.argv[4]
TF = sys.argv[5]

adata_ = sc.read_h5ad(adata_)
adata_.X = adata_.raw.X
sc.pp.normalize_total(adata_, target_sum=1e4)
sc.pp.log1p(adata_)

TF = pd.read_csv(TF, sep = '\t')
TF = TF.loc[TF['prediction'] == True, 'sequence_ID'].tolist()
adata = adata_[:, adata_.var_names.isin(TF)].copy()
adata = adata[~adata.obs[cluster_header].isin(["Mixture", 'Immune cells'])].copy()

adata.obs.rename(columns={cluster_header: cluster_header.replace(' ','_')}, inplace=True)
cluster_header = cluster_header.replace(' ','_')

ns.pp.dendrogram(adata, cluster_header, save = True, output_folder = output_folder,
                         outputfilename_suffix = label, plot = False)
adata = ns.pp.prep_medians(adata, cluster_header, use_mean = True)
adata = ns.pp.prep_binary_scores(adata, cluster_header)
# binary score stored in adata.varm['binary_scores_cell_type']

# plot exp median across cells
plt.clf()
filename = f'{output_folder}/{label}.{cluster_header}.exp_median.png'
print(f"Saving median distributions as...\n{filename}")
a = plt.figure(figsize = (6, 4))
a = plt.hist(adata.varm["medians_" + cluster_header].unstack(), bins = 100)
a = plt.title(f'{"medians_" + cluster_header} histogram')
a = plt.xlabel("medians_" + cluster_header)
a = plt.yscale("log")
a = plt.savefig(filename, bbox_inches='tight')
# plot binary score distribution across cells
plt.clf()
filename = f'{output_folder}/{label}.{cluster_header}.binary_score_distribution.png'
print(f"Saving binary_score distributions as...\n{filename}")
a = plt.figure(figsize = (6, 4))
a = plt.hist(adata.varm["binary_scores_" + cluster_header].unstack(), bins = 100)
a = plt.title(f'{"binary_scores_" + cluster_header} histogram')
a = plt.xlabel("binary_scores_" + cluster_header)
a = plt.yscale("log")
a = plt.savefig(filename, bbox_inches='tight')

results = nsforesting.NSForest(adata, cluster_header, save_supplementary = True, output_folder = output_folder,
        outputfilename_prefix = f'{label}.{cluster_header}', gene_selection = "BinaryFirst_high", 
        n_top_genes = 30, n_binary_genes = 30, n_trees = 1500)

to_plot = results.copy()
endrogram = [] # custom dendrogram order
dendrogram = list(adata.uns["dendrogram_" + cluster_header]["categories_ordered"])
to_plot["clusterName"] = to_plot["clusterName"].astype("category")
to_plot["clusterName"] = to_plot["clusterName"].cat.set_categories(dendrogram)
to_plot = to_plot.sort_values("clusterName")
to_plot = to_plot.rename(columns = {"NSForest_markers": "markers"})
to_plot.to_csv(f'{output_folder}/{label}.{cluster_header}.to_plot_results.txt', sep='\t', index=False)

markers_dict = dict(zip(to_plot["clusterName"], to_plot["markers"]))
ns.pl.dotplot(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = True, output_folder = output_folder, outputfilename_suffix = f'{label}.{cluster_header}')
ns.pl.stackedviolin(adata, markers_dict, cluster_header, dendrogram = dendrogram, save = True, output_folder = output_folder, outputfilename_suffix = f'{label}.{cluster_header}')


