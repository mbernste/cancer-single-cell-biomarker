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
    parser.add_option("-n", "--n_comps", help="Number of components")
    parser.add_option("-c", "--cluster_id", help="Cluster ID to restrict plot")
    parser.add_option("-l", "--cluster_file", help="Files storing cluster assignments. Required if '-c' flat is used.")
    (options, args) = parser.parse_args()

    dataset_metadata_f = args[0]
    out_dir = options.out_dir

    with open(dataset_metadata_f, 'r') as f:
        dataset_metadata = json.load(f)
    dataset_to_units = dataset_metadata['dataset_to_units']

    if options.n_comps:
        n_comps = int(options.n_comps)
    else:
        n_comps = 2

    sc.settings.verbosity = 3
    sc.logging.print_versions()

    for ds in sorted(set(load.DATASETS)):
        print('\nLoading data for dataset {}...'.format(ds))
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
        #sc.pp.filter_cells(ad, min_genes=10)
        sc.pp.filter_genes(ad, min_cells=3)
        if dataset_to_units[ds] == 'counts':
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad)

        print('The shape is: ', ad.X.shape)
        #print('The length of cells is: ', len(cells))
        phate_operator = phate.PHATE(n_jobs=-2, random_state=1, n_components=n_comps)
        X_phate = phate_operator.fit_transform(ad.X)

        df = pd.DataFrame(
            data=X_phate,
            index=ad.obs['cell'],
            columns=['PHATE_{}'.format(i+1) for i in range(n_comps)]
        )
        df.to_csv(
            join(out_dir, '{}_PHATE_{}.tsv'.format(ds, n_comps)),
            sep='\t'
        )


if __name__ == '__main__':
    main()

