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

    gsea_f = args[0]
    out_dir = options.out_dir

    gsea_df = pd.read_csv(gsea_f, sep='\t')
    
    plot_heatmap(gsea_df, out_dir)
    most_common_terms(gsea_df, out_dir)

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
        key=lambda x: len(term_to_tumor_clusts[x]),
        reverse=True
    )
    top_50_common_terms = common_terms_by_tumor[:50]
    print(json.dumps(
        {
            term: sorted(set([x.split('_')[0] for x in  term_to_tumor_clusts[term]]))
            for term in top_50_common_terms
        },
        indent=True
    ))

if __name__ == '__main__':
    main()

