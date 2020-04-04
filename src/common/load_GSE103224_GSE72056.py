import h5py
import numpy as np
from collections import defaultdict
import os
from os.path import join
import json

SCRIPT_DIR = os.path.dirname(__file__)
CONFIG_F = join(SCRIPT_DIR, '../../config.json')
with open(CONFIG_F, 'r') as f:
    config = json.load(f)
DATA_F = join(config['output_dir'], 'GSE103224_GSE72056.h5')
print('Data location: ', DATA_F)

# Load cell names and tumor names
with h5py.File(DATA_F, 'r') as f:
    CELLS = [
        str(x)[2:-1]
        for x in f['cell'][:]
    ]
    DATASETS = [
        str(x)[2:-1]
        for x in f['dataset'][:]
    ]
    GENE_NAMES = [
        str(x)[2:-1]
        for x in f['gene_name'][:]
    ]
    TISSUE_TYPES = [
        str(x)[2:-1]
        for x in f['tissue_type'][:]
    ]

# Map each cell to its index in the data matrix
CELL_TO_INDEX = {
    cell: index
    for index, cell in enumerate(CELLS)
}

# Map each tumor to its indices in the data matrix
DATASET_TO_INDICES = defaultdict(lambda: [])
for index, ds in enumerate(DATASETS):
    DATASET_TO_INDICES[ds].append(index)
DATASET_TO_INDICES = dict(DATASET_TO_INDICES)

def counts_matrix_all_tumors():
    with h5py.File(DATA_F, 'r') as f:
        counts = f['count'][:]
    return counts, CELLS, DATASETS

def counts_matrix_for_dataset(tumor):
    """
    Retrieve the counts matrix for a given tumor.

    Args:
        tumor: the tumor ID
    Returns:
        counts: matrix of counts for tumor
        cells: cell names correponding to the
            rows of counts
    """
    indices = DATASET_TO_INDICES[tumor]
    with h5py.File(DATA_F, 'r') as f:
        counts = f['count'][indices]
    cells = list(np.array(CELLS)[indices])
    return counts, cells


def counts_matrix_for_cells(cells):
    """
    Retrieve the counts matrix for a specific
    set of cells.

    Args:
        cells: list of N cells
    Returns:
        NxM matrix of counts for the N
        cells and M genes. Rows are returned
        in same order as input list of cells.
    """
    indices = [
        CELL_TO_INDEX[cell]
        for cell in cells
    ]
    with h5py.File(DATA_F, 'r') as f:
        counts = f['count'][indices]
    return counts


def main():
    #cells = ['PJ016_4', 'PJ016_5', 'PJ016_6']
    #print(np.sum(counts_matrix_for_cells(cells), axis=1))
    counts, X = counts_matrix_for_dataset('PJ035')
    print(list(counts[0]))
    #print(set(DATASETS))
    #print(list(np.sum(counts, axis=1)))

if __name__ == '__main__':
    main()
