import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
import h5py
import numpy as np
import json
from dash.dependencies import Input, Output
from dash_table import DataTable

from app import app
import common
import load_data
import dim_reduc
import clust_compare

@app.callback(
    Output(component_id='de-table-1', component_property='data'),
    [
        Input(component_id='select-tumor-de-1', component_property='value'),
        Input(component_id='select-cluster-de-1', component_property='value')
    ]
)
def update_de_table_1(tumor, cluster):
    genes = load_data.tumor_to_cluster_to_de_genes[tumor][cluster]
    return [
        {'de-genes-col-1': gene}
        for gene in genes
    ]


@app.callback(
    Output(component_id='de-table-2', component_property='data'),
    [   
        Input(component_id='select-tumor-de-2', component_property='value'),
        Input(component_id='select-cluster-de-2', component_property='value')
    ]
)
def update_de_table_2(tumor, cluster):
    genes = load_data.tumor_to_cluster_to_de_genes[tumor][cluster]
    return [
        {'de-genes-col-2': gene}
        for gene in genes
    ]


@app.callback(
    [
        Output(component_id='select-cluster-de-1', component_property='options'),
        Output(component_id='select-cluster-de-1', component_property='value')
    ],
    [
        Input(component_id='select-tumor-de-1', component_property='value')
    ]
)
def update_select_cluster_de_1(tumor):
    return (
        [
            {'label': 'Cluster {}'.format(clust), 'value': clust}
            for clust in load_data.tumor_to_cluster_to_de_genes[tumor].keys()
        ],
        '0'
    )


@app.callback(
    [
        Output(component_id='select-cluster-de-2', component_property='options'),
        Output(component_id='select-cluster-de-2', component_property='value')
    ],
    [
        Input(component_id='select-tumor-de-2', component_property='value')
    ]
)
def update_select_cluster_de_1(tumor):
    return (
        [
            {'label': 'Cluster {}'.format(clust), 'value': clust}
            for clust in load_data.tumor_to_cluster_to_de_genes[tumor].keys()
        ],
        '0'
    )


app.layout = html.Div(children=[
    html.H1(children='CHARTS: CHARacterizing Tumor Subpopulations'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),

    dbc.Container(fluid=True, children=[
        dcc.Tabs([
            dim_reduc.LAYOUT,
            clust_compare.LAYOUT,
            dcc.Tab(
                label='Differential Expression',
                children=[
                    dbc.Row([
                        dbc.Col([
                            html.H6("Select a tumor:"),
                            common.build_tumor_dropdown('select-tumor-de-1'),
                            html.H6("Select a cluster:"),
                            dcc.Dropdown(
                                options=[
                                    {'label': 'Cluster {}'.format(clust), 'value': clust}
                                    for clust in load_data.tumor_to_cluster_to_de_genes['PJ016'].keys()
                                ],
                                id='select-cluster-de-1',
                                value='0'
                            ),
                            DataTable(
                                columns=[{"id": "de-genes-col-1", "name": "Differentially Expressed Genes"}],
                                data=[
                                    {'de-genes-col-1': gene}
                                    for gene in load_data.tumor_to_cluster_to_de_genes['PJ016']['0']
                                ],
                                style_cell_conditional=[
                                    {
                                        'textAlign': 'center'
                                    }
                                ],
                                id='de-table-1'
                            )
                        ]),
                        dbc.Col([
                            html.H6("Select a tumor:"),
                            common.build_tumor_dropdown('select-tumor-de-2'),
                            html.H6("Select a cluster:"),
                            dcc.Dropdown(
                                options=[
                                    {'label': 'Cluster {}'.format(clust), 'value': clust}
                                    for clust in load_data.tumor_to_cluster_to_de_genes['PJ016'].keys()
                                ],
                                id='select-cluster-de-2',
                                value='0'
                            ),
                            DataTable(
                                columns=[{"id": "de-genes-col-2", "name": "Differentially Expressed Genes"}],
                                data=[
                                    {'de-genes-col-2': gene}
                                    for gene in load_data.tumor_to_cluster_to_de_genes['PJ016']['0']
                                ],
                                style_cell_conditional=[
                                    {
                                        'textAlign': 'center'
                                    }
                                ],
                                id='de-table-2'
                            )
                        ])
                    ])
                ]
            )
        ])
    ])
])

if __name__ == '__main__':
    app.run_server(debug=True)
