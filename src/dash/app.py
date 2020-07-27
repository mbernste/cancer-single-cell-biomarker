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

SCATTER_HEIGHT = '600px'

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

#app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])


############################## Load all of the data ###############################################
tumor_to_umap_df = {
    'PJ017': pd.read_csv('../../charts_data/dim_reduc/PJ017_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ016':  pd.read_csv('../../charts_data/dim_reduc/PJ016_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ018':  pd.read_csv('../../charts_data/dim_reduc/PJ018_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ035':  pd.read_csv('../../charts_data/dim_reduc/PJ035_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ025':  pd.read_csv('../../charts_data/dim_reduc/PJ025_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ030':  pd.read_csv('../../charts_data/dim_reduc/PJ030_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ048':  pd.read_csv('../../charts_data/dim_reduc/PJ048_UMAP_3.tsv', sep='\t', index_col=0),
    'PJ032':  pd.read_csv('../../charts_data/dim_reduc/PJ032_UMAP_3.tsv', sep='\t', index_col=0),
}
with open('../../charts_data/cluster/cluster_de_genes.json', 'r') as f:
    tumor_to_cluster_to_de_genes = json.load(f)

def _load_tumor_gene(tumor, gene):
    with h5py.File('GSE103224.h5', 'r') as f:
        CELLS = [
            str(x)[2:-1]
            for x in f['{}_cell'.format(tumor)][:]
        ]
        TUMORS = [
            str(x)[2:-1]
            for x in f['{}_tumor'.format(tumor)][:]
        ]
        GENE_IDS = [
            str(x)[2:-1]
            for x in f['{}_gene_id'.format(tumor)][:]
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
            expressions = np.array(f['{}_count'.format(tumor)][:,index])
        else:
            expressions = np.zeros(len(CELLS))
    df = pd.DataFrame(
        data={'expression': expressions},
        index=CELLS
    )
    return df
################################ End loading data ################################################

def _build_dim_reduc(tumor_id, gene):
    df_umap = tumor_to_umap_df[tumor_id]
    df_expr = _load_tumor_gene(tumor_id, gene)
    df = df_umap.join(df_expr)
    fig = go.Figure(data=[go.Scatter3d(
        x=df['UMAP_1'],
        y=df['UMAP_2'],
        z=df['UMAP_3'],
        mode='markers',
        marker=dict(
            size=5,
            color=df['expression'],
            colorscale='Viridis',
            opacity=0.0
        )
    )])
    return fig


@app.callback(
    Output(component_id='dim-reduc-scatter-1', component_property='figure'),
    [
        Input(component_id='select-tumor-1', component_property='value'),
        Input(component_id='color-by-gene-1', component_property='value')
    ]
)
def update_dim_reduc_1(tumor, gene):
    return _build_dim_reduc(tumor, gene)


@app.callback(
    Output(component_id='dim-reduc-scatter-2', component_property='figure'),
    [
        Input(component_id='select-tumor-2', component_property='value'),
        Input(component_id='color-by-gene-2', component_property='value')
    ]
)
def update_dim_reduc_2(tumor, gene):
    return _build_dim_reduc(tumor, gene)


@app.callback(
    Output(component_id='dim-reduc-scatter-3', component_property='figure'),
    [
        Input(component_id='select-tumor-3', component_property='value'),
        Input(component_id='color-by-gene-3', component_property='value')
    ]
)
def update_dim_reduc_3(tumor, gene):
    return _build_dim_reduc(tumor, gene)


@app.callback(
    Output(component_id='dim-reduc-scatter-4', component_property='figure'),
    [
        Input(component_id='select-tumor-4', component_property='value'),
        Input(component_id='color-by-gene-4', component_property='value')
    ]
)
def update_dim_reduc_4(tumor, gene):
    return _build_dim_reduc(tumor, gene)


@app.callback(
    Output(component_id='de-table-1', component_property='data'),
    [
        Input(component_id='select-tumor-de-1', component_property='value'),
        Input(component_id='select-cluster-de-1', component_property='value')
    ]
)
def update_de_table_1(tumor, cluster):
    genes = tumor_to_cluster_to_de_genes[tumor][cluster]
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
    genes = tumor_to_cluster_to_de_genes[tumor][cluster]
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
            for clust in tumor_to_cluster_to_de_genes[tumor].keys()
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
            for clust in tumor_to_cluster_to_de_genes[tumor].keys()
        ],
        '0'
    )


def _build_tumor_dropdown(html_id):
    return dcc.Dropdown(
        options=[
            {'label': 'PJ016 (glioma)', 'value': 'PJ016'},
            {'label': 'PJ017 (glioma)', 'value': 'PJ017'},
            {'label': 'PJ018 (glioma)', 'value': 'PJ018'},
            {'label': 'PJ025 (glioma)', 'value': 'PJ025'},
            {'label': 'PJ030 (glioma)', 'value': 'PJ030'},
            {'label': 'PJ032 (glioma)', 'value': 'PJ032'},
            {'label': 'PJ035 (glioma)', 'value': 'PJ035'},
            {'label': 'PJ048 (glioma)', 'value': 'PJ048'},
        ],
        value='PJ016',
        id=html_id
    )


app.layout = html.Div(children=[
    html.H1(children='CHARTS: CHARacterizing Tumor Subpopulations'),

    html.Div(children='''
        Dash: A web application framework for Python.
    '''),

    dbc.Container(fluid=True, children=[
        dcc.Tabs([
            dcc.Tab(
                label='Dimension Reduction',
                children=[
                    dbc.Container(fluid=True, children=[
                    dbc.Row(children=[
                        dbc.Col([
                            html.Div(children=[
                                html.H4("Plot 1"),
                                html.H6("Select a tumor to visualize:"),
                                _build_tumor_dropdown('select-tumor-1'),
                                html.H6("Color by feature (enter a gene): "),
                                dcc.Input(id='color-by-gene-1', value='OLIG1'),
                                html.H4("Plot 2"),
                                html.H6("Select a tumor to visualize:"),
                                _build_tumor_dropdown('select-tumor-2'),
                                html.H6("Color by feature (enter a gene): "),
                                dcc.Input(id='color-by-gene-2', value='OLIG1'),
                                html.H4("Plot 3"),
                                html.H6("Select a tumor to visualize:"),
                                _build_tumor_dropdown('select-tumor-3'),
                                html.H6("Color by feature (enter a gene): "),
                                dcc.Input(id='color-by-gene-3', value='OLIG1'),
                                html.H4("Plot 4"),
                                html.H6("Select a tumor to visualize:"),
                                _build_tumor_dropdown('select-tumor-4'),
                                html.H6("Color by feature (enter a gene): "),
                                dcc.Input(id='color-by-gene-4', value='OLIG1')
                            ])
                        ], width='100'),
                        dbc.Col([
                            dbc.Row([
                                dbc.Col([
                                    dbc.Row(
                                        [html.H5("Plot 1")], 
                                        align='center', 
                                        justify="center", 
                                        style={"width": "70%"}
                                    ),
                                    dbc.Row(
                                        [dcc.Graph(
                                            id='dim-reduc-scatter-1',
                                            figure=_build_dim_reduc('PJ017', 'OLIG1')
                                        )], 
                                        style={"height": SCATTER_HEIGHT, "width": "100%"}
                                    )
                                ]),
                                dbc.Col([
                                    dbc.Row(
                                        [html.H5("Plot 2")],
                                        align='center', 
                                        justify="center", 
                                        style={"width": "100%"}
                                    ), 
                                    dbc.Row(
                                        [dcc.Graph(
                                            id='dim-reduc-scatter-2',
                                            figure=_build_dim_reduc('PJ017', 'OLIG1')
                                        )],
                                        style={"height": SCATTER_HEIGHT, "width": "100%"}
                                    )
                                ])
                            ], style={"height": SCATTER_HEIGHT}),
                            dbc.Row([
                                dbc.Col([
                                    dbc.Row(
                                        [html.H5("Plot 3")],
                                        align='center',
                                        justify="center",
                                        style={"width": "70%"}
                                    ),
                                    dbc.Row(
                                        [dcc.Graph(
                                        id='dim-reduc-scatter-3',
                                        figure=_build_dim_reduc('PJ017', 'OLIG1')
                                        )],
                                        style={"height": SCATTER_HEIGHT, "width": "100%"}
                                    )
                                ]),
                                dbc.Col([
                                    dbc.Row(
                                        [html.H5("Plot 4")],
                                        align='center',
                                        justify="center",
                                        style={"width": "70%"}
                                    ),
                                    dbc.Row(
                                        [dcc.Graph(
                                        id='dim-reduc-scatter-4',
                                        figure=_build_dim_reduc('PJ017', 'OLIG1')
                                        )],
                                        style={"height": SCATTER_HEIGHT, "width": "100%"}
                                    )
                                ])
                            ], style={"height": SCATTER_HEIGHT})
                        ])
                    ])
                    ])
                ]
            ),
            dcc.Tab(
                label='Differential Expression',
                children=[
                    dbc.Row([
                        dbc.Col([
                            html.H6("Select a tumor:"),
                            _build_tumor_dropdown('select-tumor-de-1'),
                            html.H6("Select a cluster:"),
                            dcc.Dropdown(
                                options=[
                                    {'label': 'Cluster {}'.format(clust), 'value': clust}
                                    for clust in tumor_to_cluster_to_de_genes['PJ016'].keys()
                                ],
                                id='select-cluster-de-1',
                                value='0'
                            ),
                            DataTable(
                                columns=[{"id": "de-genes-col-1", "name": "Differentially Expressed Genes"}],
                                data=[
                                    {'de-genes-col-1': gene}
                                    for gene in tumor_to_cluster_to_de_genes['PJ016']['0']
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
                            _build_tumor_dropdown('select-tumor-de-2'),
                            html.H6("Select a cluster:"),
                            dcc.Dropdown(
                                options=[
                                    {'label': 'Cluster {}'.format(clust), 'value': clust}
                                    for clust in tumor_to_cluster_to_de_genes['PJ016'].keys()
                                ],
                                id='select-cluster-de-2',
                                value='0'
                            ),
                            DataTable(
                                columns=[{"id": "de-genes-col-2", "name": "Differentially Expressed Genes"}],
                                data=[
                                    {'de-genes-col-2': gene}
                                    for gene in tumor_to_cluster_to_de_genes['PJ016']['0']
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
