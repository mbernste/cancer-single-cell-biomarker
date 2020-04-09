import matplotlib as mpl
mpl.use('Agg')
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import os 
from os.path import join
import subprocess
from anndata import AnnData
import phate
from optparse import OptionParser
import h5py
import json

sys.path.append('../common')

import load_GSE103224 


def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-t", "--tumor", help="Only plot specific tumor")
    parser.add_option("-c", "--cluster_id", help="Cluster ID to restrict plot")
    parser.add_option("-n", "--n_comps", help="Number of components")
    parser.add_option("-l", "--cluster_file", help="Files storing cluster assignments. Required if '-c' flat is used.")
    (options, args) = parser.parse_args()

    in_dir = args[0]
    config_f = args[1]
    out_dir = options.out_dir

    with open(config_f, 'r') as f:
        metadata = json.load(f)
    tumors = metadata['datasets']
    tumor_to_subtype = metadata['dataset_to_subtype']
    integ_sets = metadata['integration_sets'].keys()

    if options.n_comps:
        n_comps = int(options.n_comps)
    else:
        n_comps = 2

    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    for integ_set in integ_sets:
        with h5py.File(join(in_dir, '{}_aligned.h5'.format(integ_set)), 'r') as f:
            cells = [
                str(x)[2:-1]
                for x in f['cell'][:]
            ]
            tumors = [
                str(x)[2:-1]
                for x in f['dataset'][:]
            ]
            counts = np.nan_to_num(f['count'][:])
        

        subtypes = [
            tumor_to_subtype[tumor]
            for tumor in tumors
        ]

        #cluster_df = pd.read_csv(join(in_dir, 'all_tumors_clusters_aligned.tsv'), sep='\t')
        cells_df = pd.DataFrame(data=cells, columns=['cell'], index=cells)
        #cluster_df = cells_df.join(cluster_df.set_index('cell'))
        #clusts = list(cluster_df['leiden'])
        #phate_operator = phate.PHATE(n_jobs=-2, random_state=1, n_components=n_comps)
        #X_phate = phate_operator.fit_transform(counts)

        X_phate = _run_phate(counts, n_comps)
        X_umap = _run_umap(counts, n_comps)
        _build_tsv_file(cells, tumors, subtypes, X_phate, 'PHATE', integ_set, out_dir)
        _build_tsv_file(cells, tumors, subtypes, X_umap, 'UMAP', integ_set, out_dir)

def _run_phate(counts, n_comps):
    phate_operator = phate.PHATE(
        n_jobs=-2, 
        random_state=1, 
        n_components=n_comps
    )
    return phate_operator.fit_transform(counts)
    

def _run_umap(counts, n_comps):
    ad = AnnData(
        X=counts
    )
    sc.pp.neighbors(ad, n_neighbors=15)
    sc.tl.umap(ad, n_components=n_comps)
    return ad.obsm['X_umap'] 

def _build_tsv_file(cells, tumors, subtypes, X_reduc, units, integ_set, out_dir):
        da = [
            #(cell, tumor, subtype, clust, p1, p2, p3)
            (cell, tumor, subtype, p1, p2, p3)
            #for cell, tumor, subtype, clust, (p1, p2, p3) in zip(cells, tumors, subtypes, clusts, X_phate)
            for cell, tumor, subtype, (p1, p2, p3) in zip(cells, tumors, subtypes, X_reduc)
        ]
        df = pd.DataFrame(
            data=da, 
            columns=[
                'cell', 
                'tumor', 
                'subtype', 
                #'cluster', 
                '{}1'.format(units), 
                '{}2'.format(units), 
                '{}3'.format(units)
            ]
        )
        df.to_csv(
            join(out_dir, '{}_aligned_{}_3.tsv'.format(integ_set, units)), 
            sep='\t', 
            index=False
        )

if __name__ == '__main__':
    main()

