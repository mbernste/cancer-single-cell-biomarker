#
#   Parse the counts data from the Gene Expression Omnibus for
#   accession GSE72056. Build an H5 file storing the counts.
#

from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py 
from os.path import join
from collections import defaultdict

FILE = 'GSE72056_melanoma_single_cell_revised_v2.txt'

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    root = args[0]
    out_f = args[1]

    the_genes = None
    counts = None
    tumors = []
    fname = join(root, FILE)
    print('Loading {}...'.format(fname))
    full_df = pd.read_csv(
        fname, 
        sep='\t', 
        index_col=0
    )

    raw_cells = full_df.columns 
    tumors = [
        'Mel{}'.format(int(x)).encode('utf-8') 
        for x in full_df.loc['tumor']
    ]
    df = full_df.iloc[3:]
    gene_names = df.index

    # Cast to a numpy array
    curr_counts = np.array(df).T
    print(curr_counts)
    if counts is None:
        counts = curr_counts
    else:
        counts = np.concatenate([counts, curr_counts])
    print('Current shape of the expression matrix: {}'.format(counts.shape))

    # Create more intertable cell id's
    tumor_to_count = defaultdict(lambda: 0)
    cells = []
    for tumor in tumors:
        cells.append(
            '{}_{}'.format(tumor, tumor_to_count[tumor]).encode('utf-8')
        )
        tumor_to_count[tumor] += 1

    gene_names = [
        x.encode('utf-8')
        for x in gene_names
    ]

    print('Writing to H5 file...')
    with h5py.File(out_f, 'w') as f:
        f.create_dataset('count', data=counts, compression="gzip")
        f.create_dataset('cell', data=np.array(cells))
        f.create_dataset('tumor', data=np.array(tumors))
        f.create_dataset('gene_name', data=np.array(gene_names))
    print('done.')

if __name__ == '__main__':
    main()
