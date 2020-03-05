from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py
import scanpy as sc
from anndata import AnnData

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    #parser.add_option("-a", "--a_descrip", action="store_true", help="This is a flat")
    #parser.add_option("-b", "--b_descrip", help="This is an argument")
    (options, args) = parser.parse_args()

    in_f = args[0]
    out_f = args[1]

    print('Loading data...')
    with h5py.File(in_f, 'r') as f:
        cells = [
            str(x)[2:-1].encode('utf-8')
            for x in f['cell'][:]
        ]
        tumors = [
            str(x)[2:-1].encode('utf-8')
            for x in f['tumor'][:]
        ]
        gene_ids = [
            str(x)[2:-1].encode('utf-8')
            for x in f['gene_id'][:]
        ]
        gene_names = [
            str(x)[2:-1].encode('utf-8')
            for x in f['gene_name'][:]
        ]
        counts = f['count'][:]
    print('done.')

    ad = AnnData(X=counts)
    sc.pp.normalize_total(ad, target_sum=1e6)
    sc.pp.log1p(ad)

    print('Writing to H5 file...')
    with h5py.File(out_f, 'w') as f:
        f.create_dataset('count', data=ad.X, compression="gzip")
        f.create_dataset('cell', data=np.array(cells))
        f.create_dataset('tumor', data=np.array(tumors))
        f.create_dataset('gene_id', data=np.array(gene_ids))
        f.create_dataset('gene_name', data=np.array(gene_names))
    print('done.')

if __name__ == '__main__':
    main()
