import numpy as np
import pandas as pd
import scanpy as sc
import sys
from anndata import AnnData
sys.path.append('..')

from common import load

def main():
    sc.settings.verbosity = 3
    sc.logging.print_versions()
    sc.settings.set_figure_params(dpi=80)

    counts, cells = load.counts_matrix_for_tumor('PJ016')

    assert 'MOG' in set(load.GENE_NAMES)

    ad = AnnData(
        X=counts, 
        obs=pd.DataFrame(data=cells, columns=['cell']),
        var=pd.DataFrame(
            index=load.GENE_NAMES, 
            data=load.GENE_NAMES, 
            columns=['gene_name']
        )
    )
    #print(ad_PJ048.var)

    sc.pp.normalize_total(ad, target_sum=1e6)
    sc.pp.log1p(ad)
    
    #sc.tl.pca(ad_PJ048, svd_solver='arpack')
    #sc.pl.pca(ad_PJ048, color='MOG')

    sc.pp.neighbors(ad, n_neighbors=10, n_pcs=40)
    #sc.tl.umap(ad)
    #sc.pl.umap(ad, color=['OLIG1','GFAP'])

    sc.tl.louvain(ad)
    sc.tl.paga(ad, groups='louvain')
    sc.pl.paga(ad)
    #sc.tl.draw_graph(ad, init_pos='paga')
    #sc.pl.draw_graph(ad, color=['OLIG1','GFAP'])

    sc.tl.umap(ad, init_pos='paga')
    sc.pl.umap(ad, color=['OLIG1','GFAP'])

if __name__ == '__main__':
    main()

