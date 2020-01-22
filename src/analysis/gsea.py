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
    (options, args) = parser.parse_args()

    de_genes_f = args[0]
    frac_biomarker_f = args[1]
    biomarker = args[2]
    out_dir = options.out_dir

    with open(de_genes_f, 'r') as f:
        tumor_to_cluster_to_de_genes = json.load(f)
    with open(frac_biomarker_f, 'r') as f:
        tumor_to_cluster_to_frac_nonzero = json.load(f)

    gsea_da = []
    for tumor in tumor_to_cluster_to_de_genes:
        clust_to_de_genes = tumor_to_cluster_to_de_genes[tumor]
        clust_to_frac_nonzero = tumor_to_cluster_to_frac_nonzero[tumor]

        print(clust_to_frac_nonzero)
    
        biomarker_clusts = []
        for clust in clust_to_frac_nonzero:
            de_genes = clust_to_de_genes[clust]
            if biomarker in set(de_genes) and clust_to_frac_nonzero[clust] > 0.1:
                print('{} Found in cluster {}!'.format(biomarker, clust))
                enr = gp.enrichr(
                    gene_list=de_genes,
                    gene_sets=GENE_SETS,
                    no_plot=True,
                    cutoff=0.05  # test dataset, use lower value from range(0,1)
                )
                enr.results = enr.results[enr.results["Adjusted P-value"] < GSEA_THRESH]
                print(enr.results)
                sig_terms = set(enr.results['Term'])
                gsea_da += [
                    (tumor, clust, term)
                    for term in sig_terms
                ]
    gsea_df = pd.DataFrame(
        data=gsea_da,
        columns=['tumor', 'cluster', 'GO_term']
    )
    gsea_df.to_csv(join(out_dir, 'gsea_results.tsv'), sep='\t')


if __name__ == '__main__':
    main()

