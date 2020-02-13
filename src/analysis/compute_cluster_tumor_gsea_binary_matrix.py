from optparse import OptionParser
from collections import defaultdict
import pandas as pd
import os
from os.path import join

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    (options, args) = parser.parse_args()

    gsea_f = args[0]
    out_dir = options.out_dir

    gsea_df = pd.read_csv(gsea_f, sep='\t')
    all_GO_terms = sorted(set(gsea_df['GO_term']))
    all_tumor_clusts = sorted(set(gsea_df['tumor_cluster']))

    tumor_clust_to_GO_terms = defaultdict(lambda: set())

    g = gsea_df.groupby(by='GO_term')
    mat = []
    row_terms = []
    for term, df in g:
        row = []
        term_tumor_clusts = set(df['tumor_cluster'])
        for tum_clust in all_tumor_clusts:
            row.append(int(tum_clust in term_tumor_clusts))
        mat.append(row)
        row_terms.append(term)

    binary_mat_df = pd.DataFrame(
        data=mat,
        columns=all_tumor_clusts,
        index=row_terms
    )
    binary_mat_df.to_csv(
        join(out_dir, 'gsea_binary_matrix.tsv'), 
        sep='\t'
    )

if __name__ == "__main__":
    main()
