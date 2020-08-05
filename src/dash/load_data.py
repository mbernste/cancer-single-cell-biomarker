import pandas as pd
import h5py
import numpy as np
import json

TUMOR_TO_UMAP_DF = {
    'PJ017':  pd.read_csv('../../charts_data/dim_reduc/PJ017_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ016':  pd.read_csv('../../charts_data/dim_reduc/PJ016_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ018':  pd.read_csv('../../charts_data/dim_reduc/PJ018_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ035':  pd.read_csv('../../charts_data/dim_reduc/PJ035_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ025':  pd.read_csv('../../charts_data/dim_reduc/PJ025_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ030':  pd.read_csv('../../charts_data/dim_reduc/PJ030_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ048':  pd.read_csv('../../charts_data/dim_reduc/PJ048_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ032':  pd.read_csv('../../charts_data/dim_reduc/PJ032_UMAP_3.tsv', sep='\t', index_col=0),
}

tumor_to_cell_type_df = {
    'PJ017': pd.read_csv('../../charts_data_dash/CellO_PJ017.probability.tsv', sep='\t', index_col=0),
    'PJ016': pd.read_csv('../../charts_data_dash/CellO_PJ016.probability.tsv', sep='\t', index_col=0),
    #'PJ018': pd.read_csv('../../charts_data_dash/CellO_PJ018.probability.tsv', sep='\t', index_col=0),
    #'PJ035': pd.read_csv('../../charts_data_dash/CellO_PJ035.probability.tsv', sep='\t', index_col=0),
    'PJ025': pd.read_csv('../../charts_data_dash/CellO_PJ025.probability.tsv', sep='\t', index_col=0),
    'PJ030': pd.read_csv('../../charts_data_dash/CellO_PJ030.probability.tsv', sep='\t', index_col=0),
    #'PJ048': pd.read_csv('../../charts_data_dash/CellO_PJ048.probability.tsv', sep='\t', index_col=0),
    'PJ032': pd.read_csv('../../charts_data_dash/CellO_PJ032.probability.tsv', sep='\t', index_col=0)
}

with open('../../charts_data/cluster/cluster_de_genes.json', 'r') as f:
    tumor_to_cluster_to_de_genes = json.load(f)

def load_tumor_gene_names(tumor):
    with h5py.File('../../charts.h5', 'r') as f:
        return set([
            str(x)[2:-1]
            for x in f['{}_gene_name'.format(tumor)][:]
        ])

def load_tumor_clusters(tumor):
    with h5py.File('../../charts.h5', 'r') as f:
        return set([
            str(x)
            for x in f['{}_cluster'.format(tumor)][:]
        ])

def load_tumor_clusters_for_cells(tumor):
    with h5py.File('../../charts.h5', 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tumor)][:]
        ]
        clusts = f['{}_cluster'.format(tumor)][:]
    df = pd.DataFrame(
        data={
            'cell': cells,
            'color_by': clusts
        }
    )
    df = df.set_index('cell')
    return df

def load_tumor_gene(tumor, gene):
    with h5py.File('../../charts.h5', 'r') as f:
        CELLS = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tumor)][:]
        ]
        GENE_NAMES = [
            str(x)[2:-1]
            for x in f['{}_gene_name'.format(tumor)][:]
        ]
        gene_name_to_index = {
            gene_name: index
            for index, gene_name in enumerate(GENE_NAMES)
        }
        if gene in gene_name_to_index:
            index = gene_name_to_index[gene]
            expressions = np.array(f['{}_log1_tpm'.format(tumor)][:,index])
        else:
            expressions = np.zeros(len(CELLS))
    df = pd.DataFrame(
        data={'color_by': expressions},
        index=CELLS
    )
    return df

def load_gene_mult_tumors(tum_1, tum_2, cells, gene):
    with h5py.File('../../charts.h5', 'r') as f:
        cells_1 = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tum_1)][:]
        ]
        genes_1 = [
            str(x)[2:-1]
            for x in f['{}_gene_name'.format(tum_1)][:]
        ]
        gene_name_to_index_1 = {
            gene_name: index
            for index, gene_name in enumerate(genes_1)
        }
        if gene in gene_name_to_index_1:
            index = gene_name_to_index_1[gene]
            expressions_1 = np.array(f['{}_log1_tpm'.format(tum_1)][:,index])
        else:
            expressions_1 = np.zeros(len(cells_1))
        df_1 = pd.DataFrame(
            data={'color_by': expressions_1},
            index=cells_1
        )

        cells_2 = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tum_2)][:]
        ]
        genes_2 = [
            str(x)[2:-1]
            for x in f['{}_gene_name'.format(tum_2)][:]
        ]
        gene_name_to_index_2 = {
            gene_name: index
            for index, gene_name in enumerate(genes_2)
        }
        if gene in gene_name_to_index_2:
            index = gene_name_to_index_2[gene]
            expressions_2 = np.array(f['{}_log1_tpm'.format(tum_2)][:,index])
        else:
            expressions_2 = np.zeros(len(cells_2))
        df_2 = pd.DataFrame(
            data={'color_by': expressions_2},
            index=cells_2
        ) 

        df = pd.concat([df_1, df_2])
        df = df.loc[cells]
        return df

def load_tumors_for_cells_mult_tumors(tum_1, tum_2, cells):
    tumors = []
    for cell in cells:
        if tum_1 in cell:
            tumors.append(tum_1)
        elif tum_2 in cell:
            tumors.append(tum_2)
    return pd.DataFrame(
        data={'color_by': tumors},
        index=cells
    )

def load_tumor_phate(tumor, num_dims):
    with h5py.File('../../charts.h5', 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tumor)][:]
        ]
        X_umap = f['{}_phate_{}'.format(tumor, num_dims)][:]
    df = pd.DataFrame(
        data=X_umap,
        index=cells,
        columns=['PHATE {}'.format(x+1) for x in range(num_dims)]
    )
    return df

def load_tumor_umap(tumor, num_dims):
    with h5py.File('../../charts.h5', 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tumor)][:]
        ]
        X_umap = f['{}_umap_{}'.format(tumor, num_dims)][:]
    df = pd.DataFrame(
        data=X_umap,
        index=cells,
        columns=['UMAP {}'.format(x+1) for x in range(num_dims)]
    )
    return df


