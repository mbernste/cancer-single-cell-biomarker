import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

import load_data

def build_cell_type_probability_dropdown(tumor, html_id): 
    cell_types = load_data.cell_type_probability_columns(tumor, min_prob=0.00)
    options = [
        {'label': cell_type, 'value': cell_type}
        for cell_type in sorted(cell_types)
    ]
    return dcc.Dropdown(
        options=options,
        value='endothelial cell',
        id=html_id
    )

def build_cell_type_probability_dropdown_mult_tumors(tum_1, tum_2, html_id):
    cell_types_1 = load_data.cell_type_probability_columns(tum_1, min_prob=0.00)
    cell_types_2 = load_data.cell_type_probability_columns(tum_2, min_prob=0.00)
    all_cell_types = cell_types_1 | cell_types_2
    options = [
        {'label': cell_type, 'value': cell_type}
        for cell_type in sorted(all_cell_types)
    ]
    return dcc.Dropdown(
        options=options,
        value='endothelial cell',
        id=html_id
    )

def build_hallmark_enrichment_dropdown(tumor, html_id):
    gene_sets = load_data.hallmark_gene_sets(tumor)
    options = [
        {'label': gene_set, 'value': gene_set}
        for gene_set in sorted(gene_sets)
    ]
    return dcc.Dropdown(
        options=options,
        value='HALLMARK_HYPOXIA',
        id=html_id
    )

def build_hallmark_enrichment_dropdown_mult_tumors(tum_1, tum_2, html_id):
    gene_sets_1 = set(load_data.hallmark_gene_sets(tum_1))
    gene_sets_2 = set(load_data.hallmark_gene_sets(tum_2))
    all_gene_sets = gene_sets_1 | gene_sets_2
    options = [
        {'label': gene_set, 'value': gene_set}
        for gene_set in sorted(all_gene_sets)
    ]
    return dcc.Dropdown(
        options=options,
        value='HALLMARK_HYPOXIA',
        id=html_id
    )

def build_tumor_dropdown(html_id, width=None):
    options=[
        {'label': 'PJ016 (glioma)', 'value': 'PJ016'},
        {'label': 'PJ017 (glioma)', 'value': 'PJ017'},
        {'label': 'PJ018 (glioma)', 'value': 'PJ018'},
        {'label': 'PJ025 (glioma)', 'value': 'PJ025'},
        {'label': 'PJ030 (glioma)', 'value': 'PJ030'},
        {'label': 'PJ032 (glioma)', 'value': 'PJ032'},
        {'label': 'PJ035 (glioma)', 'value': 'PJ035'},
        {'label': 'PJ048 (glioma)', 'value': 'PJ048'},
        {'label': 'LX653 (lung adenocarcinoma)', 'value': 'LX653'}, 
        {'label': 'LX661 (lung adenocarcinoma)', 'value': 'LX661'},
        {'label': 'LX675 (lung adenocarcinoma)', 'value': 'LX675'},
        {'label': 'LX676 (lung adenocarcinoma)', 'value': 'LX676'},
        {'label': 'LX679 (lung adenocarcinoma)', 'value': 'LX679'},
        {'label': 'LX680 (lung adenocarcinoma)', 'value': 'LX680'},
        {'label': 'LX682 (lung adenocarcinoma)', 'value': 'LX682'},
        {'label': 'LX684 (lung adenocarcinoma)', 'value': 'LX684'},
        {'label': 'PJ016 + PJ018', 'value': 'PJ016&PJ018'},
        {'label': 'PJ016 + PJ017', 'value': 'PJ016&PJ017'},
        {'label': 'PJ016 + PJ025', 'value': 'PJ016&PJ025'},
        {'label': 'PJ017 + PJ018', 'value': 'PJ017&PJ018'},
        {'label': 'PJ017 + PJ025', 'value': 'PJ017&PJ025'},
        {'label': 'PJ018 + PJ025', 'value': 'PJ018&PJ025'},
        {'label': 'PJ016 + LX653', 'value': 'PJ016&LX653'},
        {'label': 'PJ016 + LX676', 'value': 'PJ016&LX676'},
        {'label': 'PJ016 + LX682', 'value': 'PJ016&LX682'},
        {'label': 'PJ017 + LX653', 'value': 'PJ017&LX653'},
        {'label': 'PJ017 + LX676', 'value': 'PJ017&LX676'},
        {'label': 'PJ017 + LX682', 'value': 'PJ017&LX682'},
        {'label': 'PJ018 + LX653', 'value': 'PJ018&LX653'},
        {'label': 'PJ018 + LX676', 'value': 'PJ018&LX676'},
        {'label': 'PJ018 + LX682', 'value': 'PJ018&LX682'},
        {'label': 'PJ025 + LX653', 'value': 'PJ025&LX653'},
        {'label': 'PJ025 + LX676', 'value': 'PJ025&LX676'},
        {'label': 'PJ025 + LX682', 'value': 'PJ025&LX682'},
        {'label': 'LX653 + LX676', 'value': 'LX653&LX676'},
        {'label': 'LX653 + LX682', 'value': 'LX653&LX682'},
        {'label': 'LX676 + LX682', 'value': 'LX676&LX682'}
    ]
    if width:
        style={"width": width}
    else:
        style=None
    return dcc.Dropdown(
        options=options,
        value='PJ016',
        id=html_id,
        style=style
    )
