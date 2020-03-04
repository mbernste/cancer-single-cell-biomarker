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
    parser.add_option("-n", "--n_comps", help="Number of components")
    parser.add_option("-c", "--cluster_id", help="Cluster ID to restrict plot")
    parser.add_option("-l", "--cluster_file", help="Files storing cluster assignments. Required if '-c' flat is used.")
    (options, args) = parser.parse_args()
   
    out_dir = options.out_dir

    if options.n_comps:
        n_comps = int(options.n_comps)
    else:
        n_comps = 2

    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    for tumor in TUMORS:
        print('Loading data for tumor {}...'.format(tumor))
        counts, cells = load_GSE103224.counts_matrix_for_tumor(tumor)
        print('done.')

        ad = AnnData(
            X=counts,
            obs=pd.DataFrame(data=cells, columns=['cell']),
            var=pd.DataFrame(
                index=load_GSE103224.GENE_NAMES,
                data=load_GSE103224.GENE_NAMES,
                columns=['gene_name']
            )
        )
        sc.pp.normalize_total(ad, target_sum=1e6)
        sc.pp.log1p(ad)

        phate_operator = phate.PHATE(n_jobs=-2, random_state=1, n_components=n_comps)
        X_phate = phate_operator.fit_transform(ad.X)

        df = pd.DataFrame(
            data=X_phate,
            index=cells,
            columns=['PHATE_{}'.format(i+1) for i in range(n_comps)]
        )
        df.to_csv(
            join(out_dir, '{}_PHATE_{}.tsv'.format(tumor, n_comps)),
            sep='\t'
        )

        #np.savetxt(
        #    join(out_dir, '{}_PHATE_{}.tsv'.format(tumor, n_comps)), 
        #    X_phate, 
        #    delimiter='\t'
        #)

if __name__ == '__main__':
    main()

