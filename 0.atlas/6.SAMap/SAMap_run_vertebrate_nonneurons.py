from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
    sankey_plot, chord_plot, CellTypeTriangles,
    ParalogSubstitutions, FunctionalEnrichment,
    convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import random

random.seed(42)

fn1 = '/data/biol-brain-scSEQ/ball6395/05.SAMap_comparison/0.bin/02.atlas/Hsap.non_neurons_pr.h5ad'
fn2 = '/data/biol-brain-scSEQ/ball6395/05.SAMap_comparison/0.bin/02.atlas/Mmus.non_neurons_pr.h5ad'
fn3 = '/data/biol-brain-scSEQ/ball6395/05.SAMap_comparison/0.bin/02.atlas/Pvit.non_neurons_pr.h5ad'
fn4 = '/data/biol-brain-scSEQ/ball6395/05.SAMap_comparison/0.bin/02.atlas/Pmar.non_neurons_pr.h5ad'

sam1=SAM()
sam1.load_data(fn1)
sam2=SAM()
sam2.load_data(fn2)
sam3=SAM()
sam3.load_data(fn3)
sam4=SAM()
sam4.load_data(fn4)

sams = {'hs':sam1,'mm':sam2,'pv':sam3,'pm':sam4}
sm = SAMAP(
        sams,
        f_maps = './maps/',
        keys = {'hs':'cell_type','mm':'cell_type','pv':'cell_type','pm':'cell_type'},
        save_processed=False #if False, do not save the processed results to `*_pr.h5ad`
)


sm.run(
        pairwise=True,
        NUMITERS = 5,
        crossK = 30,
        ncpus = 8,
        hom_edge_mode = "pearson",
        N_GENE_CHUNKS = 6
        )
#sm.run(neigh_from_keys = {'hs':True,'mm':True,'pv':True}) #
#samap = sm.samap # SAM object with three species stitched togetherfilenamesfilenames

import pickle
with open('vertebrate_non_neurons.pkl', 'wb') as f:  # open a text file
    pickle.dump(sm, f)
f.close()
