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

FRAC_THRESH = 0.1

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
    in_dir = args[1]
    tumor_to_clust_to_de_genes_f = args[2]
    out_dir = options.out_dir

    # Load data
    with open(tumor_to_clust_to_de_genes_f, 'r') as f:
        tumor_to_clust_to_de_genes = json.load(f)

    # Compute the fraction of 
    tumor_to_clust_to_frac_expr = _fraction_of_cells_expressing_biomarker_per_cluster(
        biomarker, 
        in_dir  
    )

    # Determine which clusters to keep
    tumor_clusters = []
    for tumor in TUMORS:
        for clust, frac_nonzero in tumor_to_clust_to_frac_expr[tumor].items(): 
            if frac_nonzero > FRAC_THRESH \
                and biomarker in tumor_to_clust_to_de_genes[tumor][str(clust)]:
                    tumor_clusters.append((tumor, str(clust)))

    # Write output
    with open(join(out_dir, 'cluster_fraction_biomarker.json'), 'w') as f:
        json.dump(
            tumor_to_clust_to_frac_expr,
            f, 
            indent=4
        )
    df_kept_clusts = pd.DataFrame(
        data=tumor_clusters,
        columns=['tumor', 'cluster']
    )
    df_kept_clusts.to_csv(
        join(out_dir, 'clusters_expressing_biomarker.tsv'), 
        sep='\t'
    )

def _fraction_of_cells_expressing_biomarker_per_cluster(
        biomarker, 
        in_dir 
    ):
    tumor_to_clust_to_frac_expr = defaultdict(lambda: {})
    for tumor in TUMORS:
        print('Loading data for tumor {}...'.format(tumor))
        counts, cells = load_GSE103224.counts_matrix_for_tumor(tumor)
        print('done.')

        cluster_df = pd.read_csv(
            join(in_dir, '{}_clusters.tsv'.format(tumor)),
            sep='\t'
        )
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
        biomarker_ind = load_GSE103224.GENE_NAMES.index(biomarker)
        grouped = grouped = cluster_df.groupby('louvain') 
        for clust, group in grouped:
            indices = [int(x) for x in group.index]
            clust_X = ad.X[indices]
            biomarker_X = clust_X[:,biomarker_ind]
            frac_nonzero = sum([1 for x in biomarker_X if x > 0]) / len(biomarker_X)
            tumor_to_clust_to_frac_expr[tumor][clust] = frac_nonzero
    return tumor_to_clust_to_frac_expr

def plot_heatmap(gsea_df, out_dir):
    all_sig_terms = sorted(set(gsea_df['GO_term']))
    tumor_clust_to_terms = defaultdict(lambda: set())
    for index, row in gsea_df.iterrows():
        tumor = row['tumor']
        clust = row['cluster']
        term = row['GO_term']
        tumor_clust = '{}_{}'.format(tumor, clust)
        tumor_clust_to_terms[tumor_clust].add(term)

    tumor_clusts = sorted(tumor_clust_to_terms.keys())
    heatmap_da = [
        [
            int(term in tumor_clust_to_terms[tumor_clust])
            for tumor_clust in tumor_clusts
        ]
        for term in all_sig_terms
    ]        
    heatmap_df = pd.DataFrame(
        data=heatmap_da, 
        columns=tumor_clusts, 
        index=all_sig_terms
    )

    res = sns.clustermap(heatmap_df, xticklabels=True, yticklabels=False, row_cluster=True, cmap='Greys')
    res.savefig(
        join(out_dir, 'common_GO_terms_heatmap.png'),
        format='png',
        dpi=150
    )

def most_common_terms(gsea_df, out_dir):
    term_to_tumor_clusts = defaultdict(lambda: set())
    for index, row in gsea_df.iterrows():
        tumor = row['tumor']
        clust = row['cluster']
        term = row['GO_term']
        tumor_clust = '{}_{}'.format(tumor, clust)
        term_to_tumor_clusts[term].add(tumor_clust)

    common_terms = sorted(
        term_to_tumor_clusts.keys(),
        key=lambda x: len(term_to_tumor_clusts[x]),
        reverse=True
    )
    top_50_common_terms = common_terms[:50]
    #print(json.dumps(
    #    {
    #        term: sorted(term_to_tumor_clusts[term])
    #        for term in top_50_common_terms
    #    }, 
    #    indent=True
    #))

    common_terms_by_tumor = sorted(
        term_to_tumor_clusts.keys(),
        key=lambda x: len(set([y.split('_')[0] for y in term_to_tumor_clusts[x]])),
        reverse=True
    )
    top_50_common_terms = common_terms_by_tumor[:50]

    with open(join(out_dir, 'most_common_GO_terms.json'), 'w') as f: 
        json.dump(
            {
                term: sorted(term_to_tumor_clusts[term])
                for term in top_50_common_terms
            },
            f,
            indent=True
        )

if __name__ == '__main__':
    main()

