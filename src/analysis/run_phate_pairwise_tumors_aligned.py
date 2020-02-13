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

sys.path.append('../common')

import load_GSE103224 


TUMORS = [
    'PJ016',
    'PJ018',
    'PJ048',
    'PJ030',
    'PJ035',
    'PJ025',
    'PJ017',
    'PJ032'
]

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
    out_dir = options.out_dir

    if options.n_comps:
        n_comps = int(options.n_comps)
    else:
        n_comps = 2

    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    for t_i, tumor_i in enumerate(TUMORS):
        for t_j, tumor_j in enumerate(TUMORS):
            if t_j == t_i:
                continue
            #try:
            fname = '{}_{}.h5'.format(tumor_i, tumor_j)
            print('Loading {}'.format(fname))        
            with h5py.File(join(in_dir, fname), 'r') as f:
                print(f.keys())
                cells = [
                    str(x)[2:-1]
                    for x in f['cell'][:]
                ]
                tumors = [
                    str(x)[2:-1]
                    for x in f['tumor'][:]
                ]
                counts = f['count'][:]       
            #print(cells)
            #except: # TODO REMOVE THIS
            #    continue

            phate_operator = phate.PHATE(n_jobs=-2, random_state=1, n_components=n_comps)
            X_phate = phate_operator.fit_transform(counts)

            da = [
                (cell, tumor, p1, p2, p3)
                for cell, tumor, (p1, p2, p3) in zip(cells, tumors, X_phate)
            ]
            df = pd.DataFrame(data=da, columns=['cell', 'tumor', 'PHATE1', 'PHATE2', 'PHATE3'])
            df.to_csv(join(out_dir, '{}_{}_aligned_PHATE_3.tsv'.format(tumor_i, tumor_j)), sep='\t', index=False)

if __name__ == '__main__':
    main()

