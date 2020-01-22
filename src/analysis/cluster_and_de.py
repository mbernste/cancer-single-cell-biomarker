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

    biomarker = args[0]
    out_dir = options.out_dir

    if options.resolution:
        resolution = float(options.resolution)
    else:
        resolution = 0.8
    
    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    tumor_to_biomarker_clusts = defaultdict(lambda: [])
    tumor_to_clust_to_de_genes = defaultdict(lambda: {})
    tumor_to_clust_to_sig_terms = defaultdict(lambda: {})
    tumor_to_clust_to_frac_biomarker = defaultdict(lambda: {})
    all_sig_terms = set()

    gsea_da = []
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
        biomarker_ind = load_GSE103224.GENE_NAMES.index(biomarker)
        #print(ad_PJ048.var)

        sc.pp.normalize_total(ad, target_sum=1e6)
        sc.pp.log1p(ad)
        sc.pp.pca(ad, n_comps=50)
        sc.pp.neighbors(ad)
        sc.tl.louvain(ad, resolution=resolution)

        # Record cluster of each cell
        ad.obs[['cell', 'louvain']].to_csv(
            join(out_dir, '{}_clusters.tsv'.format(tumor)), 
            sep='\t',
            index=False
        )

        grouped = ad.obs.groupby('louvain')
        clust_to_frac_nonzero = {}
        for clust, group in grouped:
            indices = [int(x) for x in group.index]
            clust_X = ad.X[indices]
            biomarker_X = clust_X[:,biomarker_ind]
            print(clust)
            print(len(biomarker_X))
            frac_nonzero = sum([1 for x in biomarker_X if x > 0]) / len(biomarker_X)
            print(frac_nonzero)
            clust_to_frac_nonzero[clust] = frac_nonzero
            tumor_to_clust_to_frac_biomarker[tumor][clust] = frac_nonzero

        sc.tl.rank_genes_groups(
            ad, groupby='louvain', method='wilcoxon', n_genes=len(load_GSE103224.GENE_NAMES)
        )
        
        biomarker_clusts = []
        for clust in sorted(set(ad.obs['louvain'])):
            gene_to_pval = {
                gene: pval
                for gene, pval in zip(
                    ad.uns['rank_genes_groups']['names'][clust], 
                    ad.uns['rank_genes_groups']['pvals_adj'][clust]
                )
            }
            de_genes = [
                gene
                for gene, pval in gene_to_pval.items()
                if pval < 0.05
            ]
            tumor_to_clust_to_de_genes[tumor][clust] = de_genes
    with open(join(out_dir, 'cluster_de_genes.json'), 'w') as f:
        json.dump(tumor_to_clust_to_de_genes, f, indent=True)
    with open(join(out_dir, 'cluster_fraction_biomarker.json'), 'w') as f:
        json.dump(tumor_to_clust_to_frac_biomarker, f, indent=True)

        """
            if biomarker in gene_to_pval \
                and gene_to_pval[biomarker] < 0.05 \
                and clust_to_frac_nonzero[clust] > 0.1:
                print('{} Found in cluster {}!'.format(biomarker, clust))
                biomarker_clusts.append(clust)
                tumor_to_biomarker_clusts[tumor].append(clust)

                enr = gp.enrichr(
                    gene_list=de_genes,
                    gene_sets=GENE_SETS,
                    no_plot=True,
                    cutoff=0.05  # test dataset, use lower value from range(0,1)
                )
                enr.results = enr.results[enr.results["Adjusted P-value"] < GSEA_THRESH]
                sig_terms = set(enr.results['Term'])
                tumor_to_clust_to_sig_terms[clust][clust] = sig_terms
                all_sig_terms.update(sig_terms)
                gsea_da += [
                    (tumor, clust, term)
                    for term in sig_terms
                ]
        """
    #with open(join(out_dir, 'biomarker_clusters.json'), 'w') as f:
    #    json.dump(tumor_to_biomarker_clusts, indent=True)

    gsea_df = pd.DataFrame(
        data=gsea_da,
        columns=['tumor', 'cluster', 'GO_term']
    )
    gsea_df.to_csv(join(out_dir, 'gsea_results.tsv'), sep='\t')
    
    

    #da = [
    #    ()
    #]
    #for sig_term in all_sig_terms:

    #print('Tumor to clusters expressing biomarker: {}'.format(tumor_to_biomarker_clusts))
    #print('Tumor to cluster to sig terms: {}'.format(tumor_to_clust_to_sig_terms))
    #print(ad.uns['rank_genes_groups']['names']['0'])
    #sc.pl.rank_genes_groups(ad)

def _create_plots(ad, out_dir, genes):
    for gene in genes:
        # Cell type markers
        fig, ax = plt.subplots(1,1,figsize=(6,6))
        ax = sc.pl.scatter(ad, x='PHATE 1', y='PHATE 2', color=gene, ax=ax, legend_loc='none', show=False)
        ax.set_xticks([])
        ax.set_yticks([])
        l, b, w, h = fig.axes[-1].get_position().bounds
        ll, bb, ww, hh = fig.axes[0].get_position().bounds
        fig.axes[0].set_position([ll-0.05, bb, ww, hh])
        fig.axes[-1].set_position([ll+ww-0.05, b, w, h])
        #plt.tight_layout()
        fig.savefig(
            join(out_dir, '{}.png'.format(gene)),
            format='png',
            dpi=150
            #bbox_inches='tight'
        )


if __name__ == '__main__':
    main()

