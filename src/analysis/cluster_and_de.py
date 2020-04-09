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
import json

sys.path.append('..')

#from common import load_GSE103224_GSE72056 as load 
from common import load_GSE103224_GSE123904 as load


def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-r", "--resolution", help="Resolution")
    (options, args) = parser.parse_args()

    dataset_metadata_f = args[0]
    with open(dataset_metadata_f, 'r') as f:
        dataset_metadata = json.load(f)
    dataset_to_units = dataset_metadata['dataset_to_units']

    out_dir = options.out_dir

    if options.resolution:
        resolution = float(options.resolution)
    else:
        resolution = 0.8
    
    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    ds_to_clust_to_de_genes = defaultdict(lambda: {})

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
        sc.pp.pca(ad, n_comps=50)
        sc.pp.neighbors(ad)
        sc.tl.leiden(ad, resolution=resolution)

        # Record cluster of each cell
        res_str = str(resolution).replace('.', '_')
        ad.obs = ad.obs.rename(columns={'leiden': 'cluster'})
        out_f = join(
            out_dir, 
            '{}_clusters.res_{}.tsv'.format(
                ds, 
                res_str
            )
        )
        print('Writing data to file ', out_f)
        ad.obs[['cell', 'cluster']].to_csv(
            out_f,
            sep='\t',
            index=False
        )

        # Rank genes within each cluster
        sc.tl.rank_genes_groups(
            ad, 
            groupby='cluster', 
            method='wilcoxon', 
            n_genes=len(load.GENE_NAMES)
        )
        for clust in sorted(set(ad.obs['cluster'])):
            gene_to_pval = {
                gene: pval
                for gene, pval in zip(
                    ad.uns['rank_genes_groups']['names'][clust], 
                    ad.uns['rank_genes_groups']['pvals_adj'][clust]
                )
            }
            gene_to_logfold = {
                gene: lf
                for gene, lf in zip(
                    ad.uns['rank_genes_groups']['names'][clust],
                    ad.uns['rank_genes_groups']['logfoldchanges'][clust]
                )
            }
            de_genes = [
                gene
                for gene in gene_to_pval.keys()
                if gene_to_pval[gene] < 0.01 
                and gene_to_logfold[gene] > 1
            ]
            ds_to_clust_to_de_genes[ds][clust] = de_genes
    with open(join(out_dir, 'cluster_de_genes.json'), 'w') as f:
        json.dump(ds_to_clust_to_de_genes, f, indent=True)


if __name__ == '__main__':
    main()

