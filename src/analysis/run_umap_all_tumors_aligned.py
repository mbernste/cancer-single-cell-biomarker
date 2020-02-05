import matplotlib as mpl
#mpl.use('Agg')
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

sys.path.append('../common')

import load_GSE103224 


TUMORS = [
    'PJ016',
    'PJ018',
    'PJ048',
    'PJ030',
    'PJ025',
    'PJ035',
    'PJ017',
    'PJ032'
]

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-t", "--tumor", help="Only plot specific tumor")
    parser.add_option("-c", "--cluster_id", help="Cluster ID to restrict plot")
    parser.add_option("-l", "--cluster_file", help="Files storing cluster assignments. Required if '-c' flat is used.")
    (options, args) = parser.parse_args()

    in_dir = args[0]
    out_dir = options.out_dir

    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    with h5py.File(join(in_dir, 'GSE103224_aligned.h5'), 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['cell'][:]
        ]
        tumors = [
            str(x)[2:-1]
            for x in f['tumor'][:]
        ]
        genes = [
            str(x)[2:-1]
            for x in f['gene_name'][:]
        ]
        counts = f['count'][:]

    ad = AnnData(
        X=counts,
        obs=pd.DataFrame(
            data=[x for x in zip(cells, tumors)], 
            columns=['cell', 'tumor']
        ),
        var=pd.DataFrame(
            index=genes,
            data=genes,
            columns=['gene_name']
        )
    )
    sc.tl.pca(ad, n_comps=50)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)

    da = [
        (cell, tumor, u1, u2) 
        for cell, tumor, (u1, u2) in zip(cells, tumors, ad.obsm['X_umap'])
    ]
    df = pd.DataFrame(data=da, columns=['cell', 'tumor', 'UMAP1', 'UMAP2'])
    #sc.pl.umap(ad, color='tumor')
    #plt.show()

    df.to_csv(join(out_dir, 'all_tumors_aligned_UMAP.tsv'), sep='\t', index=False)

if __name__ == '__main__':
    main()

