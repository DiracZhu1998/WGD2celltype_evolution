# import dependencies
import sys, os, glob, re, pickle
from functools import partial
from collections import OrderedDict
import operator as op
from cytoolz import compose

import pandas as pd
import seaborn as sns
import numpy as np
import scanpy as sc
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt

from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss

from IPython.display import HTML, display

from save_logos import plotLogo

# Set maximum number of jobs for Scanpy.
sc.settings.njobs = 4
ncores = 10

MAIN = sys.argv[1] + '/'
RESOURCES_FOLDERNAME = "/mnt/data01/yuanzhen/01.Vertebrate_cell_evo/09.GRNs_detection/1.mouse/0.bin/"
RESULTS_FOLDERNAME = MAIN + "/1.results/"
FIGURES_FOLDERNAME = MAIN + "/2.figures/"
sc.settings.figdir = FIGURES_FOLDERNAME
BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"
PYSCENIC_PATH = "/home/zoo/ball6395/software/miniconda3/envs/scenicpy/bin/pyscenic"
os.system(f'mkdir {RESULTS_FOLDERNAME} {FIGURES_FOLDERNAME}')

TFS_FNAME = os.path.join(RESOURCES_FOLDERNAME, 'allTFs_mm.txt')
RANKING_DBS_FNAMES = list(map(lambda fn: os.path.join(RESOURCES_FOLDERNAME, "cisTargets/", fn),
    ['mm9-500bp-upstream-7species.mc9nr.feather','mm9-tss-centered-10kb-7species.mc9nr.feather']))
print(RANKING_DBS_FNAMES)
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDERNAME, 'motifs-v9-nr.mgi-m0.001-o0.0.tbl')
CELL_ANNOTATIONS_FNAME = os.path.join(MAIN, "0.cellInfo.subset.csv")
EXP_MTX_COUNTS_FNAME = os.path.join(MAIN, '0.exprMat.subset.csv')

# Results file name
DATASET_ID = "Mouse_brain"
METADATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.metadata.csv'.format(DATASET_ID))
EXP_MTX_QC_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.qc.log_umi.csv'.format(DATASET_ID))
ADJACENCIES_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.adjacencies.tsv'.format(DATASET_ID))
MOTIFS_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.modules.csv'.format(DATASET_ID))
MOTIFS_LOGO = os.path.join(RESULTS_FOLDERNAME, '{}.module_motifs_logo.html'.format(DATASET_ID))
REGULONS_DAT_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.regulons.dat'.format(DATASET_ID))
REGULONS_CSV_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.regulons_df.csv'.format(DATASET_ID))
AUCELL_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.auc.csv'.format(DATASET_ID))
AUCELL_RSS_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.auc_rss.celltype_family.csv'.format(DATASET_ID))
AUCELL_RSS_C_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.auc_rss.iterative_cluster.csv'.format(DATASET_ID))
BIN_MTX_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.bin.csv'.format(DATASET_ID))
THR_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.thresholds.csv'.format(DATASET_ID))
ANNDATA_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.h5ad'.format(DATASET_ID))
LOOM_FNAME = os.path.join(RESULTS_FOLDERNAME, '{}.loom'.format(DATASET_ID))

# STEP 0: Load data and preprocessing
df_metadata = pd.read_csv(CELL_ANNOTATIONS_FNAME)
df_umi = pd.read_csv(EXP_MTX_COUNTS_FNAME, index_col=0)
print(df_umi.shape)
# filter df_umi based on geneFiltering function R version
nCountsPerGene = df_umi.sum(axis=1, skipna=True)
nCellsPerGene = (df_umi > 0).sum(axis=1, skipna=True)
minCountsPerGene = 1 * 0.01 * df_umi.shape[1] # equal to at least 1% cells express that gene with 3UMI or 3% cells express that gene with 1 UMI 
minSamples = df_umi.shape[1] * 0.005 # at least 1 % cells express
genesLeft_minReads = nCountsPerGene[nCountsPerGene > minCountsPerGene].index
genesLeft_minCells = nCellsPerGene[nCellsPerGene > minSamples].index
df_umi = df_umi.loc[genesLeft_minReads.intersection(genesLeft_minCells)]
print(df_umi.shape)

adata = sc.AnnData(X=df_umi.T.sort_index())
df_obs = df_metadata[['CellID','DonorID','TaxonomyRank3', 'TaxonomyRank4', 'Species', 'Refined family', 'Refined subtype']].set_index('CellID').sort_index()
adata.obs = df_obs
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Store non-log transformed data as raw. This data can be used via the use_raw parameters available for many functions.
# In the scanpy's tutorials this is used to stored all genes in log-transformed counts before retaining only Highly Variable Genes (HVG). 
# Because in this case no filtering is done we use this feature to store raw counts.
adata.raw = adata
sc.pp.log1p(adata)
adata
adata.write_h5ad(ANNDATA_FNAME)

df_log_umi = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
df_log_umi.to_csv(EXP_MTX_QC_FNAME)

# STEP 1: Network inference based on GRNBoost2, calculate co-expression adjacencies
# Output: List of adjacencies between a TF and its targets stored in ADJACENCIES_FNAME
os.system(f'{PYSCENIC_PATH} grn {EXP_MTX_QC_FNAME} {TFS_FNAME} -o  {ADJACENCIES_FNAME} --num_workers {ncores}')

# STEP 2-3: Construct module and prune module to regulon
DBS_PARAM = ' '.join(RANKING_DBS_FNAMES)
os.system(f'echo {PYSCENIC_PATH} ctx  {ADJACENCIES_FNAME} {DBS_PARAM}  --annotations_fname {MOTIF_ANNOTATIONS_FNAME}  --expression_mtx_fname {EXP_MTX_QC_FNAME}  --output {MOTIFS_FNAME}  --num_workers {ncores}')
os.system(f'{PYSCENIC_PATH} ctx  {ADJACENCIES_FNAME} {DBS_PARAM}  --annotations_fname {MOTIF_ANNOTATIONS_FNAME}  --expression_mtx_fname {EXP_MTX_QC_FNAME}  --output {MOTIFS_FNAME}  --num_workers {ncores}')

df_motifs = load_motifs(MOTIFS_FNAME)

# plot logo and motif
plotLogo(df = df_motifs,  output = MOTIFS_LOGO, base_url = BASE_URL)

regulons = df2regulons(df_motifs)
with open(REGULONS_DAT_FNAME, "wb") as f:
    pickle.dump(regulons, f)

regulons_df = pd.concat(
            [pd.DataFrame.from_dict(reg.gene2weight, orient="index", columns=["value"]).assign(TF=reg.name)
                     for reg in regulons], ignore_index=False)
regulons_df.to_csv(REGULONS_CSV_FNAME)

#STEP 4: Cellular enrichment aka AUCell
auc_mtx = aucell(df_log_umi, regulons, num_workers = ncores)
auc_mtx.to_csv(AUCELL_MTX_FNAME)

auc_mtx = pd.read_csv(AUCELL_MTX_FNAME)
auc_mtx.set_index("Cell", inplace=True )
rss_cellType = regulon_specificity_scores(auc_mtx, df_metadata.loc[auc_mtx.index]['Refined family'])
rss_cellType.T.to_csv(AUCELL_RSS_FNAME)
rss_cluster = regulon_specificity_scores(auc_mtx, df_metadata.loc[auc_mtx.index]['Refined subtype'])
rss_cluster.T.to_csv(AUCELL_RSS_C_FNAME)



# STEP 5 - Regulon activity binarization
bin_mtx, thresholds = binarize(auc_mtx, num_workers=20)
bin_mtx.to_csv(BIN_MTX_FNAME)
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv(THR_FNAME)


