import subprocess
import sys
from anndata import AnnData
import scanpy as sc
import phate
import h5py
import pandas as pd
import numpy as np
from optparse import OptionParser

def main():
    usage = "" # TODO
    parser = OptionParser(usage=usage)
    parser.add_option("-o", "--out_dir", help="Directory to write output")
    parser.add_option("-w", "--overwrite", action="store_true", help="Overwrite")
    (options, args) = parser.parse_args()

    tum_1 = args[0]
    tum_2 = args[1]
    h5_f = args[2]
    overwrite = options.overwrite

    if not overwrite:
        with h5py.File(h5_f, 'r') as f:
            umap_2_key = '{}&{}_umap_2'.format(tum_1, tum_2)
            umap_3_key = '{}&{}_umap_3'.format(tum_1, tum_2)
            phate_2_key = '{}&{}_phate_2'.format(tum_1, tum_2)
            phate_3_key = '{}&{}_phate_3'.format(tum_1, tum_2)
            cells_key = '{}&{}_cell'.format(tum_1, tum_2)
            curr_keys = set(f.keys())
            targ_keys = set([
                umap_2_key,
                umap_3_key,
                phate_2_key,
                phate_3_key,
                cells_key,
            ])
            if len(targ_keys - curr_keys) == 0:
                print("The dataset was detected in the HDF5 file. Not processing this data.")
                return

    subprocess.run('mkdir tmp', shell=True)
    cmd = 'Rscript integration.R {f} {t1} {t2} ./tmp/aligned.{t1}_{t2}.h5'.format(f=h5_f, t1=tum_1, t2=tum_2)
    subprocess.run(cmd, shell=True)

    with h5py.File('./tmp/aligned.{}_{}.h5'.format(tum_1, tum_2), 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['{}.{}_cell'.format(tum_1, tum_2)][:]
        ]
        expression = f['{}.{}_log1_tpm'.format(tum_1, tum_2)][:]
    print(cells)
    ad = AnnData(
        X=expression,
        obs=pd.DataFrame(
            data={'cell': cells}, 
            index=cells
        )
    )

    X_phate_2 = run_phate(ad, 2)
    X_umap_2 = run_umap(ad, 2)
    X_phate_3 = run_phate(ad, 3)
    X_umap_3 = run_umap(ad, 3)

    with h5py.File(h5_f, 'r+') as f:
        try:
            del f['{}&{}_cell'.format(tum_1, tum_2)]
        except KeyError:
            pass
        f.create_dataset(
            '{}&{}_cell'.format(tum_1, tum_2),
            data=np.array([x.encode('utf-8') for x in cells]),
            compression="gzip"
        )
        try:
            del f['{}&{}_umap_2'.format(tum_1, tum_2)]
        except KeyError:
            pass
        f.create_dataset(
            '{}&{}_umap_2'.format(tum_1, tum_2),
            data=np.array(X_umap_2, dtype=np.float32),
            compression="gzip"
        )
        try:
            del f['{}&{}_phate_2'.format(tum_1, tum_2)]
        except KeyError:
            pass
        f.create_dataset(
            '{}&{}_phate_2'.format(tum_1, tum_2),
            data=np.array(X_phate_2, dtype=np.float32),
            compression="gzip"
        )
        try:
            del f['{}&{}_umap_3'.format(tum_1, tum_2)]
        except KeyError:
            pass
        f.create_dataset(
            '{}&{}_umap_3'.format(tum_1, tum_2),
            data=np.array(X_umap_3, dtype=np.float32),
            compression="gzip"
        )
        try:
            del f['{}&{}_phate_3'.format(tum_1, tum_2)]
        except KeyError:
            pass
        f.create_dataset(
            '{}&{}_phate_3'.format(tum_1, tum_2),
            data=np.array(X_phate_3, dtype=np.float32),
            compression="gzip"
        )

def run_phate(ad, n_comps):
    print("Running PHATE with {} dimensions...".format(n_comps))
    phate_operator = phate.PHATE(
        n_jobs=-2,
        random_state=1,
        n_components=n_comps
    )
    X_phate = phate_operator.fit_transform(ad.X)
    return X_phate


def run_umap(ad, n_comps):
    print("Running UMAP with {} dimensions...".format(n_comps))
    sc.tl.pca(ad, n_comps=50)
    if ad.X.shape[0] < 300:
        n_neighbors = int(ad.X.shape[0] * 0.05)
    else:
        n_neighbors = 15
    sc.pp.neighbors(ad, n_neighbors=n_neighbors)
    sc.tl.umap(ad, n_components=n_comps)
    return ad.obsm['X_umap']

if __name__ == '__main__':
    main()
