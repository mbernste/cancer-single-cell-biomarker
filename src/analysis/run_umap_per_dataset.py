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
import json

sys.path.append('..')

from common import load_GSE103224_GSE72056 as load

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-t", "--tumor", help="Only plot specific tumor")
    parser.add_option("-n", "--n_comps", help="Number of components")
    parser.add_option("-l", "--cluster_file", help="Files storing cluster assignments. Required if '-c' flat is used.")
    (options, args) = parser.parse_args()

    if options.n_comps:
        n_comps = int(options.n_comps)
    else:
        n_comps = 2

    dataset_metadata_f = args[0]
    out_dir = options.out_dir

    with open(dataset_metadata_f, 'r') as f:
        dataset_metadata = json.load(f)
    dataset_to_units = dataset_metadata['dataset_to_units']

    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    for ds in sorted(set(load.DATASETS)):
        print('Loading data for dataset {}...'.format(ds))
        counts, cells = load.counts_matrix_for_dataset(ds)
        print('done.')

        ad = AnnData(
            X=counts,
            obs=pd.DataFrame(data=cells, columns=['cell']),
            var=pd.DataFrame(
                index=load.GENE_NAMES,
                data=load.GENE_NAMES,
                columns=['gene_name']
            )
        )
        if dataset_to_units[ds] == 'counts':
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad)
        sc.tl.pca(ad, n_comps=50)
        if ad.X.shape[0] < 300:
            n_neighbors = int(ad.X.shape[0] * 0.05)
        else:
            n_neighbors = 15
        sc.pp.neighbors(ad, n_neighbors=n_neighbors)
        sc.tl.umap(ad, n_components=n_comps)

        df = pd.DataFrame(
            data=ad.obsm['X_umap'],
            index=cells,
            columns=['UMAP_{}'.format(i+1) for i in range(n_comps)]
        )
        df.to_csv(
            join(out_dir, '{}_UMAP_{}.tsv'.format(ds, n_comps)),
            sep='\t'
        )

if __name__ == '__main__':
    main()

