from samalg import SAM
import pandas as pd
import scanpy as sc

import sys, re

data = sys.argv[1]
k = int(sys.argv[2])
npcs = int(sys.argv[3])
n_genes = int(sys.argv[4]) # like 5000
label = sys.argv[5] # like Hsap_neurons
label = f'{label}.{k}.{npcs}.{n_genes}.h5ad'


sam=SAM()
sam.load_data(data)

sam.preprocess_data()
sam.run(stopping_condition=5e-3,
        k = k, npcs = npcs,
        n_genes = n_genes,
        seed = 42)

sam.save_anndata(label)

