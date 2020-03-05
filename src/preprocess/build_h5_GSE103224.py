#
#   Parse the counts data from the Gene Expression Omnibus for
#   accession GSE103224. Build an H5 file storing the counts.
#

from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py 
from os.path import join

TUMOR_TO_FILE = {
    'PJ016': 'GSM2758471_PJ016.filtered.matrix.txt',
    'PJ018': 'GSM2758473_PJ018.filtered.matrix.txt',
    'PJ030': 'GSM2758475_PJ030.filtered.matrix.txt',
    'PJ035': 'GSM2758477_PJ035.filtered.matrix.txt',
    'PJ017': 'GSM2758472_PJ017.filtered.matrix.txt',
    'PJ025': 'GSM2758474_PJ025.filtered.matrix.txt',
    'PJ032': 'GSM2758476_PJ032.filtered.matrix.txt',
    'PJ048': 'GSM2940098_PJ048.filtered.matrix.txt'
}

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    root = args[0]
    out_f = args[1]

    print('Finding genes common to all datasets...')
    the_genes = None
    counts = None
    cells = []
    tumors = []
    for tumor, tumor_fname in TUMOR_TO_FILE.items():
        tumor_f = join(root, tumor_fname)
        print('Loading {}...'.format(tumor_f))
        df = pd.read_csv(
            tumor_f, 
            sep='\t', 
            header=None, 
            index_col=0
        )

        # Extract the genes
        genes = list(df.index)
        if the_genes is None:
            the_genes = list(genes)
            gene_names = list(df.iloc[:,0])
        assert frozenset(genes) == frozenset(the_genes)

        # Drop the gene-names column
        keep_cols = df.columns[1:]
        df = df[keep_cols]

        # Re-order rows according to the global gene 
        # ordering
        df = df.loc[the_genes]

        # Cast to a numpy array
        curr_counts = np.array(df).T
        if counts is None:
            counts = curr_counts
        else:
            counts = np.concatenate([counts, curr_counts])
        print('Current shape of the expression matrix: {}'.format(counts.shape))

        cells += [
            '{}_{}'.format(tumor, i).encode('utf-8')
            for i in np.arange(curr_counts.shape[0])
        ]
        tumors += [
            tumor.encode('utf-8')
            for i in np.arange(curr_counts.shape[0])
        ]


    the_genes = [
        x.encode('utf-8')
        for x in the_genes
    ]
    gene_names = [
        x.encode('utf-8')
        for x in gene_names
    ]
    print('done.')

    print('Writing to H5 file...')
    with h5py.File(out_f, 'w') as f:
        f.create_dataset('count', data=counts, compression="gzip")
        f.create_dataset('cell', data=np.array(cells))
        f.create_dataset('tumor', data=np.array(tumors))
        f.create_dataset('gene_id', data=np.array(the_genes))
        f.create_dataset('gene_name', data=np.array(gene_names))
    print('done.')

if __name__ == '__main__':
    main()
