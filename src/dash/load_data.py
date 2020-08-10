import pandas as pd
import h5py
import numpy as np
import json

with open('../../charts_data/cluster/cluster_de_genes.json', 'r') as f:
    tumor_to_cluster_to_de_genes = json.load(f)

def load_tumor_gene_names(tumor):
    with h5py.File('../../charts.h5', 'r') as f:
        return set([
            str(x)[2:-1]
            for x in f['{}_gene_name'.format(tumor)][:]
        ])

def cell_type_probability_columns(tumor, min_prob=None):
    with h5py.File('../../charts.h5', 'r') as f:
        cols = [
            str(x)[2:-1]
            for x in f['{}_cell_type_probability_columns'.format(tumor)]
        ]
        probs = f['{}_cell_type_probability'.format(tumor)][:]
        mins = np.min(probs, axis=0) 
        df = pd.DataFrame(
            data={'probability': mins},
            index=cols
        )
        return set(df.loc[df['probability'] > min_prob].index)

def hallmark_gene_sets(tumor):
    with h5py.File('../../charts.h5', 'r') as f:
        cols = [
            str(x)[2:-1]
            for x in f['{}_hallmark_gene_set_name'.format(tumor)]
        ]
    return sorted(set(cols))

def load_tumor_cell_type_classifications(tumor):
    with h5py.File('../../charts.h5', 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tumor)][:]
        ]
        cell_types = [
            str(x)[2:-1]
            for x in f['{}_predicted_cell_type'.format(tumor)]
        ]
    df = pd.DataFrame(
        data={
            'cell': cells,
            'color_by': cell_types
        }
    )
    df = df.set_index('cell')
    return df

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
        cells = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tumor)][:]
        ]
        gene_names = [
            str(x)[2:-1]
            for x in f['{}_gene_name'.format(tumor)][:]
        ]
        gene_name_to_index = {
            gene_name: index
            for index, gene_name in enumerate(gene_names)
        }
        if gene in gene_name_to_index:
            index = gene_name_to_index[gene]
            expressions = np.array(f['{}_log1_tpm'.format(tumor)][:,index])
        else:
            expressions = np.zeros(len(cells))
    df = pd.DataFrame(
        data={'color_by': expressions},
        index=cells
    )
    return df


def load_tumor_cell_type_probabilities(tumor, cell_type):
    with h5py.File('../../charts.h5', 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tumor)][:]
        ]
        cell_types = [
            str(x)[2:-1]
            for x in f['{}_cell_type_probability_columns'.format(tumor)][:]
        ]
        cell_type_to_index = {
            cell_type: index
            for index, cell_type in enumerate(cell_types)
        }
        index = cell_type_to_index[cell_type]
        probabilities = np.array(f['{}_cell_type_probability'.format(tumor)][:,index])
    df = pd.DataFrame(
        data={'color_by': probabilities},
        index=cells
    )
    return df


def load_tumor_hallmark_enrichment(tumor, gene_set):
    with h5py.File('../../charts.h5', 'r') as f:
        cells = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tumor)][:]
        ]
        gene_sets = [
            str(x)[2:-1]
            for x in f['{}_hallmark_gene_set_name'.format(tumor)][:]
        ]
        gene_set_to_index = {
            gene_set: index
            for index, gene_set in enumerate(gene_sets)
        }
        index = gene_set_to_index[gene_set]
        scores = np.array(f['{}_hallmark_gsva'.format(tumor)][:,index])
    df = pd.DataFrame(
        data={'color_by': scores},
        index=cells
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


def load_tumor_cell_type_probabilities_mult_tumors(
        tum_1,
        tum_2,
        cells,
        cell_type
    ):
    with h5py.File('../../charts.h5', 'r') as f:
        cells_1 = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tum_1)][:]
        ]
        cell_types_1 = [
            str(x)[2:-1]
            for x in f['{}_cell_type_probability_columns'.format(tum_1)][:]
        ]
        cell_type_to_index_1 = {
            cell_type: index
            for index, cell_type in enumerate(cell_types_1)
        }
        if cell_type in cell_type_to_index_1:
            index = cell_type_to_index_1[cell_type]
            probs_1 = np.array(f['{}_cell_type_probability'.format(tum_1)][:,index])
        else:
            probs_1 = np.zeros(len(cells_1))
        df_1 = pd.DataFrame(
            data={'color_by': probs_1},
            index=cells_1
        )

        cells_2 = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tum_2)][:]
        ]
        cell_types_2 = [
            str(x)[2:-1]
            for x in f['{}_cell_type_probability_columns'.format(tum_2)][:]
        ]
        cell_type_to_index_2 = {
            cell_type: index
            for index, cell_type in enumerate(cell_types_2)
        }
        if cell_type in cell_type_to_index_2:
            index = cell_type_to_index_2[cell_type]
            probs_2 = np.array(f['{}_cell_type_probability'.format(tum_2)][:,index])
        else:
            probs_2 = np.zeros(len(cells_2))
        df_2 = pd.DataFrame(
            data={'color_by': probs_2},
            index=cells_2
        )

        df = pd.concat([df_1, df_2])
        df = df.loc[cells]
        return df


def load_cell_type_classifications_mult_tumors(
        tum_1,
        tum_2,
        cells
    ):
    with h5py.File('../../charts.h5', 'r') as f:
        cells_1 = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tum_1)][:]
        ]
        cell_types_1 = [
            str(x)[2:-1]
            for x in f['{}_predicted_cell_type'.format(tum_1)][:]
        ]
        df_1 = pd.DataFrame(
            data={'color_by': cell_types_1},
            index=cells_1
        )

        cells_2 = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tum_2)][:]
        ]
        cell_types_2 = [
            str(x)[2:-1]
            for x in f['{}_predicted_cell_type'.format(tum_2)][:]
        ]
        df_2 = pd.DataFrame(
            data={'color_by': cell_types_2},
            index=cells_2
        )

        df = pd.concat([df_1, df_2])
        df = df.loc[cells]
        return df

def load_clusters_mult_tumors(
        tum_1,
        tum_2,
        cells
    ):
    with h5py.File('../../charts.h5', 'r') as f:
        cells_1 = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tum_1)][:]
        ]
        clusters_1 = [
            '{}_{}'.format(tum_1, x)
            for x in f['{}_cluster'.format(tum_1)][:]
        ]
        df_1 = pd.DataFrame(
            data={'color_by': clusters_1},
            index=cells_1
        )

        cells_2 = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tum_2)][:]
        ]
        clusters_2 = [
            '{}_{}'.format(tum_2, x)
            for x in f['{}_cluster'.format(tum_2)][:]
        ]
        df_2 = pd.DataFrame(
            data={'color_by': clusters_2},
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


