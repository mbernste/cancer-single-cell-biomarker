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

from common import load_GSE103224_GSE72056 as load 

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-r", "--resolution", help="Resolution")
    (options, args) = parser.parse_args()

    in_dir = args[0]
    tumor_to_clust_to_de_genes_f = args[1]
    out_dir = options.out_dir

    # Load the tumor-metadata
    with open(join(in_dir, 'tumor_metadata.json'), 'r') as f:
        metadata = json.load(f)
    dataset_to_units = metadata['dataset_to_units']

    # Load data
    tumor_to_cluster_df = {
        ds: pd.read_csv(
            join(in_dir, 'cluster', '{}_clusters.res_0_8.tsv'.format(ds)),
            sep='\t'
        )
        for ds in sorted(set(load.DATASETS))
    }
    with open(tumor_to_clust_to_de_genes_f, 'r') as f:
        tumor_to_clust_to_de_genes = json.load(f)

    # Load the tumor-clusters
    the_tumor_clusters = []
    for ds in sorted(set(load.DATASETS)):
        for cluster in set(tumor_to_cluster_df[ds]['cluster']):
            the_tumor_clusters.append((ds, cluster))
    the_tumor_clusters = sorted(the_tumor_clusters)

    # Determine which genes to consider
    the_genes = set()
    for ds, clust_to_de_genes in tumor_to_clust_to_de_genes.items():
        for clust, de_genes in clust_to_de_genes.items():
            the_genes.update(de_genes)
    the_genes = sorted(the_genes)
    print("{} total DE genes".format(len(the_genes)))
    gene_to_index = {
        gene: ind
        for ind, gene in enumerate(load.GENE_NAMES)
    }

    # Compute the fraction of
    tumor_to_ad = {}
    tumor_to_cluster_to_inds = defaultdict(lambda: {})
    da = []
    for ds in sorted(set(load.DATASETS)):
        counts, cells = load.counts_matrix_for_dataset(ds)
        ad = AnnData(
            X=counts,
            obs=pd.DataFrame(
                data=cells,
                columns=['cell']
            ),
            var=pd.DataFrame(
                index=load.GENE_NAMES,
                data=load.GENE_NAMES,
                columns=['gene_name']
            )
        )
        if dataset_to_units[ds] == 'counts':
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad)

        grouped = tumor_to_cluster_df[ds].groupby('cluster')
        for clust, group in grouped:
            print("Examining tumor {}, cluster {}".format(ds, clust))
            row = ["{}_{}".format(ds, clust)]
            indices = [int(x) for x in group.index]
            X_clust = clust_X = ad.X[indices]
            for gene in the_genes:
                X_gene = X_clust.T[gene_to_index[gene]]
                frac_nonzero = len([x for x in X_gene if x > 0]) / len(X_gene)
                if gene in tumor_to_clust_to_de_genes[ds][str(clust)]:
                    row.append(frac_nonzero)
                else:
                    row.append(0.0)
            da.append(row)
    df = pd.DataFrame(data=da, columns=['tumor_cluster']+the_genes)
    df.to_csv(join(out_dir, 'tumor_cluster_gene_expression.tsv'), sep='\t', index=False)


def _fraction_of_cells_expressing_biomarker_per_cluster(
        gene,
        gene_to_index,
        in_dir,
        tumor_to_ad,
        tumor_to_cluster_to_inds
    ):
    print('Computing fractions for gene {}'.format(gene))
    tumor_to_clust_to_frac_expr = defaultdict(lambda: {})
    for tumor in sorted(set(load.DATASETS)):
        ad = tumor_to_ad[tumor]
        for clust, indices in tumor_to_cluster_to_inds[tumor].items():
            clust_X = ad.X[indices]
            gene_X = clust_X[:,gene_to_index[gene]]
            frac_nonzero = len([x for x in gene_X if x > 0]) / len(gene_X)
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

