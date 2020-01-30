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
    (options, args) = parser.parse_args()
   
    in_dir = args[0]
    dim_reduc = args[1]
    out_dir = options.out_dir

    cells = load_GSE103224.CELLS
    tumors = load_GSE103224.TUMORS

    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    phate_X = np.loadtxt(
        join(in_dir, 'all_tumors_{}.tsv'.format(dim_reduc)),
        delimiter='\t'
    )
    print(phate_X.shape)
    #phate_X = phate_X.T

    ad = AnnData(
        X=phate_X,
        obs=pd.DataFrame(
            data=[x for x in zip(cells, tumors)], 
            columns=['cell', 'tumor']
        ),
        var=pd.DataFrame(
            data=[
                '{} 1'.format(dim_reduc), 
                '{} 2'.format(dim_reduc)
            ],
            index=[
                '{} 1'.format(dim_reduc), 
                '{} 2'.format(dim_reduc)
            ],
            columns=['dimension']
        )
    )

    print(ad)

    # Color points by cluster
    fig, ax = plt.subplots(1,1,figsize=(8,6))
    ax = sc.pl.scatter(
        ad, 
        x='{} 1'.format(dim_reduc), 
        y='{} 2'.format(dim_reduc), 
        color='tumor', 
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
        join(out_dir, 'all_tumors_{}.png'.format(dim_reduc)),
        format='png',
        dpi=150
    )

if __name__ == '__main__':
    main()

