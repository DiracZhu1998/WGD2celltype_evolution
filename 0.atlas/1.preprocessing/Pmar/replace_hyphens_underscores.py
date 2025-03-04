import anndata as ad

import sys, re

file_ = sys.argv[1]
h5ad = ad.read_h5ad(file_)

gene_do_not_replace = []

gene_do_not_replace = []
new_gene_names = [
        gene if gene in gene_do_not_replace else gene.replace('-', '_')
        for gene in h5ad.var.index
]

h5ad.var.index = new_gene_names

h5ad.write_h5ad(file_,  compression='gzip')
