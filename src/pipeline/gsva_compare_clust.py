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
from collections import defaultdict
import gseapy as gp
import h5py
import run_gsva

sys.path.append('..')

GENE_SETS = ['GO_Biological_Process_2018']
#GENE_SETS = ['GO_Molecular_Function_2018']

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-w", "--overwrite", action='store_true', help="Overwrite data in the HDF5 file if there's a dataset already present")
    (options, args) = parser.parse_args()

    h5_f = args[0]
    overwrite = options.overwrite
    out_dir = options.out_dir

    the_tumors = set()
    with h5py.File(h5_f, 'r') as f:
        for k in f.keys():
            the_tumors.add(k.split('_')[0])
        # The symbol '&' indicates it is an integrated dataset
        the_tumors = sorted([x for x in the_tumors if '&' not in x])

    print(the_tumors)

    gene_set_to_genes = run_gsva._parse_gene_sets('./gene_sets/h.all.v7.1.symbols.gmt')

    # Get the columns from the first tumor and use as a reference
    # ensuring that all of tumor's have the same column order
    with h5py.File(h5_f, 'r') as f:
        the_cols = [
            str(x)[2:-1]
            for x in f['{}_hallmark_gene_set_name'.format(list(the_tumors)[0])]
        ]
    print('{} total columns'.format(len(the_cols)))


    the_clusts = []
    mat = []
    for tumor in the_tumors:
        print("Retrieving cluster scores for tumor {}".format(tumor))
        with h5py.File(h5_f, 'r') as f:
            cells = [
                str(x)[2:-1]
                for x in f['{}_cell'.format(tumor)][:]
            ]
            cols = [
                str(x)[2:-1]
                for x in f['{}_hallmark_gene_set_name'.format(tumor)]
            ]
            scores = f['{}_hallmark_gsva'.format(tumor)][:]
            clusts = f['{}_cluster'.format(tumor)][:]

        # Map each cluster to its cells
        clust_to_cells = defaultdict(lambda: [])
        for cell, clust in zip(cells, clusts):
            clust_to_cells[clust].append(cell)
        clust_to_cells = dict(clust_to_cells)

        df = pd.DataFrame(
            data=scores,
            columns=cols,
            index=cells
        )
        df = df[the_cols]
        
        for clust, cells in clust_to_cells.items():
            scores = np.array(df.loc[cells])
            means = np.mean(scores, axis=0)
            the_clusts.append('{}_{}'.format(tumor, clust))
            mat.append(means)
    mat = np.array(mat)

    print(the_clusts)
    print(mat.shape)

    with h5py.File(h5_f, 'r+') as f:
        try:
            del f['gsva_compare_cluster']
        except KeyError:
            pass
        try:
            del f['gsva_compare_cluster_gene_set_name']
        except KeyError:
            pass
        try:
            del f['gsva_compare_cluster_cluster']
        except KeyError:
            pass
        f.create_dataset(
            'gsva_compare_cluster',
            data=mat,
            compression='gzip'
        )
        f.create_dataset(
            'gsva_compare_cluster_gene_set_name',
            data=np.array([
                x.encode('utf-8')
                for x in the_cols
            ]),
            compression='gzip'
        )
        f.create_dataset(
            'gsva_compare_cluster_cluster',
            data=np.array([
                x.encode('utf-8')
                for x in the_clusts
            ]),
            compression='gzip'
        )

if __name__ == '__main__':
    main()

