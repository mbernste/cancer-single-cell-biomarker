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

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-t", "--tumor", help="Only plot specific tumor")
    parser.add_option("-c", "--cluster_id", help="Cluster ID to restrict plot")
    parser.add_option("-l", "--cluster_file", help="Files storing cluster assignments. Required if '-c' flat is used.")
    (options, args) = parser.parse_args()
   
    gene = args[0]
    in_dir = args[1]
    dim_reduc = args[2]
    out_dir = options.out_dir

    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    assert gene in set(load_GSE103224.GENE_NAMES)

    X_dr = np.loadtxt(
        join(in_dir, 'all_tumors_{}.tsv'.format(dim_reduc)),
        delimiter='\t'
    )
    
    counts, cells, tumors = load_GSE103224.counts_matrix_all_tumors()

    ad = AnnData(
        X=counts,
        obs=pd.DataFrame(
            data=[
                (cell, tumor, x1, x2)
                for cell, tumor, (x1, x2) in zip(cells, tumors, X_dr)
            ],
            columns=[
                'cell', 
                'tumor', 
                '{} 1'.format(dim_reduc), 
                '{} 2'.format(dim_reduc)
            ]
        ),
        var=pd.DataFrame(
            data=load_GSE103224.GENE_NAMES,
            index=load_GSE103224.GENE_NAMES,
            columns=['gene']
        )
    )

    # Color points by cluster
    fig, ax = plt.subplots(1,1,figsize=(8,6))
    ax = sc.pl.scatter(
        ad, 
        x='{} 1'.format(dim_reduc), 
        y='{} 2'.format(dim_reduc), 
        color=gene, 
        ax=ax, 
        legend_loc='right margin', 
        show=False
    )
    ax.set_xticks([])
    ax.set_yticks([])
    l, b, w, h = fig.axes[-1].get_position().bounds
    ll, bb, ww, hh = fig.axes[0].get_position().bounds
    plt.tight_layout()
    fig.savefig(
        join(out_dir, 'all_tumors_{}_{}.png'.format(dim_reduc, gene)),
        format='png',
        dpi=200
        #bbox_inches='tight'
    )

if __name__ == '__main__':
    main()

