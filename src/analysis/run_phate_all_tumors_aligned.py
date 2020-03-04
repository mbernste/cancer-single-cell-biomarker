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
        config = json.load(f)
    tumors = config['tumors']
    tumor_to_subtype = config['tumor_to_subtype']

    if options.n_comps:
        n_comps = int(options.n_comps)
    else:
        n_comps = 2

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

    subtypes = [
        tumor_to_subtype[tumor]
        for tumor in tumors
    ]

    cluster_df = pd.read_csv(join(in_dir, 'all_tumors_clusters_aligned.tsv'), sep='\t')
    cells_df = pd.DataFrame(data=cells, columns=['cell'], index=cells)
    cluster_df = cells_df.join(cluster_df.set_index('cell'))
    clusts = list(cluster_df['leiden'])

    #ad = AnnData(
    #    X=X,
    #    obs=pd.DataFrame(data=[x for x in zip(cells, tumors)], columns=['cell', 'tumor']),
    #    var=pd.DataFrame(
    #        index=genes,
    #        data=genes,
    #        columns=['gene_name']
    #    )
    #)
    phate_operator = phate.PHATE(n_jobs=-2, random_state=1, n_components=n_comps)
    X_phate = phate_operator.fit_transform(counts)

    #np.savetxt(
    #    join(out_dir, 'all_tumors_aligned_PHATE_{}.tsv'.format(n_comps)), 
    #    X_phate, 
    #    delimiter='\t'
    #)
    #ad.obs.to_csv(join(out_dir, 'all_tumors_aligned_PHATE_annot_{}.tsv'.format(n_comps)), sep='\t', index=False)

    da = [
        (cell, tumor, subtype, clust, p1, p2, p3)
        for cell, tumor, subtype, clust, (p1, p2, p3) in zip(cells, tumors, subtypes, clusts, X_phate)
    ]
    df = pd.DataFrame(data=da, columns=['cell', 'tumor', 'subtype', 'cluster', 'PHATE1', 'PHATE2', 'PHATE3'])
    #sc.pl.umap(ad, color='tumor')
    #plt.show()

    df.to_csv(join(out_dir, 'all_tumors_aligned_PHATE_3.tsv'), sep='\t', index=False)

if __name__ == '__main__':
    main()

