import sys
from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py
import scanpy as sc
from anndata import AnnData

sys.path.append('../common')

#import load_GSE103224_GSE72056 as load
import load_GSE103224_GSE123904 as load

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    #parser.add_option("-a", "--a_descrip", action="store_true", help="This is a flat")
    #parser.add_option("-b", "--b_descrip", help="This is an argument")
    (options, args) = parser.parse_args()

    metadata_f = args[0]
    out_f = args[1]

    with open(metadata_f, 'r') as f:
        metadata = json.load(f)
    dataset_to_units = metadata['dataset_to_units']
    
    counts, cells, datasts = load.counts_matrix_all_datasets()
    print('Number of genes before filtering: ', len(load.GENE_NAMES))
    print('Loading full dataset...')
    ad = AnnData(
        X=counts,
        obs=pd.DataFrame(data=cells, columns=['cell']),
        var=pd.DataFrame(
            index=load.GENE_NAMES,
            data=load.GENE_NAMES,
            columns=['gene_name']
        )
    )
    print('done.')
    sc.pp.filter_genes(ad, min_cells=30)
    the_genes = ad.var.index
    print('Number of genes after filtering: ', len(the_genes))

    all_counts = None
    all_cells = []
    all_datasets = []
    for ds in sorted(set(load.DATASETS)):
        if ds not in metadata['datasets']:
            continue
        print('\nLoading data for dataset {}...'.format(ds))
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
        ad = ad[:,the_genes]
        if dataset_to_units[ds] == 'counts':
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad)
        if all_counts is None:
            all_counts = ad.X
        else:
            print(all_counts.shape)
            print(ad.X.shape)
            all_counts = np.concatenate([all_counts, ad.X])
        datasets = [
            ds
            for i in range(len(cells))
        ]
        all_datasets += datasets
        all_cells += cells

    print('The full matrix has shape ', all_counts.shape)
    all_cells = _encode_txt(all_cells)
    all_datasets = _encode_txt(all_datasets)
    the_genes = _encode_txt(the_genes)

    print('Writing to H5 file...')
    with h5py.File(out_f, 'w') as f:
        f.create_dataset('count', data=all_counts, compression="gzip")
        f.create_dataset('cell', data=np.array(all_cells))
        f.create_dataset('dataset', data=np.array(all_datasets))
        f.create_dataset('gene_name', data=np.array(the_genes))
    print('done.')

def _encode_txt(l):
    return [x.encode('utf-8') for x in l]

if __name__ == '__main__':
    main()
