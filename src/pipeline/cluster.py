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
import h5py

sys.path.append('..')

GENE_SETS = ['GO_Biological_Process_2018']
#GENE_SETS = ['GO_Molecular_Function_2018']

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-r", "--resolution", help="Resolution")
    parser.add_option("-w", "--overwrite", action="store_true", help="Overwrite existing clustering data in the database")
    (options, args) = parser.parse_args()

    h5_f = args[0]
    out_dir = options.out_dir
    overwrite = options.overwrite

    if options.resolution:
        resolution = float(options.resolution)
    else:
        resolution = 1.0
    
    sc.settings.verbosity = 3
    sc.logging.print_versions()
    
    the_tumors = set()
    with h5py.File(h5_f, 'r') as f:
        for k in f.keys():
            the_tumors.add(k.split('_')[0])
        # The symbol '&' indicates it is an integrated dataset
        the_tumors = sorted([x for x in the_tumors if '&' not in x])

    print(the_tumors)

    for tumor in the_tumors:
        # Determine whether to cluster this tumor or not
        if not overwrite:
            with h5py.File(h5_f, 'r') as f:
                if '{}_cluster'.format(tumor) in f.keys():
                    print('Found tumor {} in the database. Skipping clustering.'.format(tumor))
                    continue

        print("Clustering tumor {}".format(tumor))
        with h5py.File(h5_f, 'r') as f:
            cells = [
                str(x)[2:-1]
                for x in f['{}_cell'.format(tumor)][:]
            ]
            expression = f['{}_log1_tpm'.format(tumor)][:]
        ad = AnnData(
            X=expression, 
            obs=pd.DataFrame(data=cells, columns=['cell'], index=cells)
        )
        sc.pp.pca(ad, n_comps=50)
        sc.pp.neighbors(ad)
        sc.tl.leiden(ad, resolution=resolution)

        clusters = ad.obs['leiden']
        assert tuple(ad.obs.index) == tuple(cells)
        
        with h5py.File(h5_f, 'r+') as f:
            try:
                del f['{}_cluster'.format(tumor)]
            except KeyError:
                pass
            f.create_dataset(
                '{}_cluster'.format(tumor), 
                data=np.array(clusters, dtype=np.int32), 
                compression="gzip"
            )
if __name__ == '__main__':
    main()

