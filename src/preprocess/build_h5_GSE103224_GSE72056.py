from optparse import OptionParser
import pandas as pd
import json
import numpy as np
import h5py 
import sys
from anndata import AnnData
import scanpy as sc
import magic

sys.path.append('..')

from common import load_GSE103224
from common import load_GSE72056

#ROOT = '/Users/matthewbernstein/Development/single-cell-hackathon/data/GSE103224_RAW'
OUT_F = '/Users/matthewbernstein/Development/single-cell-hackathon/data/GSE103224_GSE72056.h5'

def main():
    usage = "" # TODO 
    parser = OptionParser(usage=usage)
    #parser.add_option("-a", "--a_descrip", action="store_true", help="This is a flat")
    #parser.add_option("-b", "--b_descrip", help="This is an argument")
    (options, args) = parser.parse_args()

    with open('gene_synonyms.json', 'r') as f:
        syn_to_sym = json.load(f)

    print('Loading GSE72056 data...')
    GSE72056_counts, GSE72056_cells, GSE72056_tumors = load_GSE72056.counts_matrix_all_tumors()
    print('done.')
    #print('Loading GSE103224 data.')
    run_magic = False
    if not run_magic:
        GSE103224_counts, GSE103224_cells, GSE103224_tumors = load_GSE103224.counts_matrix_all_tumors()
    else:
        GSE103224_counts = None
        GSE103224_tumors = []
        GSE103224_cells = []
        #for ds in sorted(set(load_GSE103224.TUMORS)):
        for ds in ['PJ035']:
            print('Loading tumor {} from dataset GSE103224...'.format(ds))
            counts, cells = load_GSE103224.counts_matrix_for_tumor(ds)
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
            sc.pp.normalize_total(ad, target_sum=1e6)
            sc.pp.log1p(ad)
            GSE103224_cells += cells
            magic_operator = magic.MAGIC()
            X = magic_operator.fit_transform(ad.X)
            print(list(X[0]))
            if GSE103224_counts is None:
                GSE103224_counts = X
            else:
                GSE103224_counts = np.concatenate([GSE103224_counts, X])
            GSE103224_tumors += [
                ds
                for i in np.arange(len(cells))
            ]

    # Load gene information for GSE72056
    GSE72056_genes = []
    for gene in load_GSE72056.GENE_NAMES:
        if gene in syn_to_sym:
            GSE72056_genes.append(syn_to_sym[gene])
        else:
            GSE72056_genes.append(gene)
    GSE72056_gene_to_index = {
        gene: index
        for index, gene in enumerate(GSE72056_genes)
    }
    
    # Create tissue types for GSE72056
    GSE72056_tissue_types = [
        'melanoma'
        for x in GSE72056_tumors
    ]

    # Load gene information for GSE103224
    GSE103224_genes = []
    for gene in load_GSE103224.GENE_NAMES:
        if gene in syn_to_sym:
            GSE103224_genes.append(syn_to_sym[gene])
        else:
            GSE103224_genes.append(gene)
    GSE103224_gene_to_index = {
        gene: index
        for index, gene in enumerate(GSE103224_genes)
    }
    GSE103224_tissue_types = [
        'glioma'
        for x in GSE103224_tumors
    ]

    print('{} genes in GSE103224.'.format(len(GSE103224_genes)))
    print('{} genes in GSE72056.'.format(len(GSE72056_genes)))

    common_genes = sorted(set(GSE103224_genes) & set(GSE72056_genes))
    print('{} genes in common.'.format(len(common_genes)))
    print('{} genes in GSE72056 not in GSE103224:\n{}'.format(
        len(set(GSE72056_genes) - set(GSE103224_genes)), 
        set(GSE72056_genes) - set(GSE103224_genes)
    ))

    GSE72056_gene_inds = [
        GSE72056_gene_to_index[gene]
        for gene in common_genes
    ]
    GSE103224_gene_inds = [
        GSE103224_gene_to_index[gene]
        for gene in common_genes
    ]

    GSE72056_counts = GSE72056_counts[:,GSE72056_gene_inds]
    GSE103224_counts = GSE103224_counts[:,GSE103224_gene_inds]

    #print(list(GSE103224_counts[0]))
    print(list(np.sum(GSE103224_counts, axis=1)))

    print('Shape of GSE103224 counts matrix: ', GSE103224_counts.shape)
    print('Shape of GSE72056 counts matrix: ', GSE72056_counts.shape)
    counts = np.concatenate([GSE103224_counts, GSE72056_counts])
    cells = GSE103224_cells + GSE72056_cells
    data_sets = GSE103224_tumors + GSE72056_tumors
    tissue_types = GSE103224_tissue_types + GSE72056_tissue_types
    print('Shape of full counts matrix: ', counts.shape)

    common_genes = [
        x.encode('utf-8')
        for x in common_genes
    ]
    cells = [
        x.encode('utf-8')
        for x in cells
    ]
    data_sets = [
        x.encode('utf-8')
        for x in data_sets
    ]
    tissue_types = [
        x.encode('utf-8')
        for x in tissue_types
    ]

    print('Writing to H5 file...')
    with h5py.File(OUT_F, 'w') as f:
        f.create_dataset('count', data=counts, compression="gzip")
        f.create_dataset('cell', data=np.array(cells))
        f.create_dataset('dataset', data=np.array(data_sets))
        f.create_dataset('gene_name', data=np.array(common_genes))
        f.create_dataset('tissue_type', data=np.array(tissue_types))
    print('done.')


if __name__ == '__main__':
    main()
