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

from common import load_GSE103224 

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

#TUMORS = [
#    'PJ030',
#    'PJ035'
#]

GENE_SETS = ['GO_Biological_Process_2018']
#GENE_SETS = ['GO_Molecular_Function_2018']

GSEA_THRESH = 0.05

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-r", "--resolution", help="Resolution")
    (options, args) = parser.parse_args()

    out_dir = options.out_dir

    if options.resolution:
        resolution = float(options.resolution)
    else:
        resolution = 0.8
    
    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    tumor_to_clust_to_de_genes = defaultdict(lambda: {})

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
        #ad.raw = ad
        sc.pp.log1p(ad)
        sc.pp.pca(ad, n_comps=50)
        sc.pp.neighbors(ad)
        sc.tl.louvain(ad, resolution=resolution)

        # Record cluster of each cell
        res_str = str(resolution).replace('.', '_')
        ad.obs[['cell', 'louvain']].to_csv(
            join(out_dir, '{}_clusters.res_{}.tsv'.format(
                tumor, 
                res_str
            )), 
            sep='\t',
            index=False
        )

        # Rank genes within each cluster
        sc.tl.rank_genes_groups(
            ad, 
            groupby='louvain', 
            method='wilcoxon', 
            n_genes=len(load_GSE103224.GENE_NAMES)
        )
        for clust in sorted(set(ad.obs['louvain'])):
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
            tumor_to_clust_to_de_genes[tumor][clust] = de_genes
    with open(join(out_dir, 'cluster_de_genes.json'), 'w') as f:
        json.dump(tumor_to_clust_to_de_genes, f, indent=True)


if __name__ == '__main__':
    main()

