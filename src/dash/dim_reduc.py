import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc
import plotly.express as px
import plotly.graph_objects as go
from dash.dependencies import Input, Output

from app import app
import load_data
import common

SCATTER_HEIGHT = '550px'

FIG_DIM = 500

STYLE = 'legendpoints'

# Color blind palette from:
# https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
PALETTE = [
    "#004949",
    "#009292",
    "#ff6db6",
    "#ffb6db",
    "#490092",
    "#006ddb",
    "#b66dff",
    "#6db6ff",
    "#b6dbff",
    "#920000",
    "#924900",
    "#db6d00",
    "#24ff24",
    "#ffff6d",
    "#000000"
]


@app.callback(
    Output(component_id='dim-reduc-scatter-1', component_property='figure'),
    [   
        Input(component_id='select-tumor-1', component_property='value'),
        Input(component_id='dim-reduc-alg-1', component_property='value'),
        Input(component_id='color-by-feature-1', component_property='value'),
        Input(component_id='select-feature-category-1', component_property='value')
    ]
)
def update_dim_reduc_1(tumor, algo, gene, category):
    return _build_dim_reduc(tumor, algo, gene, category)


@app.callback(
    Output(component_id='color-by-feature-container-1', component_property='children'),
    [
        Input(component_id='select-tumor-1', component_property='value'),
        Input(component_id='select-feature-category-1', component_property='value')
    ]
)
def update_feature_category_selector_1(tumor, category):
    return build_features_selector('color-by-feature-1', tumor, category)

@app.callback(
    Output(component_id='dim-reduc-scatter-2', component_property='figure'),
    [
        Input(component_id='select-tumor-2', component_property='value'),
        Input(component_id='dim-reduc-alg-2', component_property='value'),
        Input(component_id='color-by-feature-2', component_property='value')
    ]
)
def update_dim_reduc_2(tumor, algo, gene):
    return _build_dim_reduc(tumor, algo, gene, 'gene')


@app.callback(
    Output(component_id='dim-reduc-scatter-3', component_property='figure'),
    [
        Input(component_id='select-tumor-3', component_property='value'),
        Input(component_id='dim-reduc-alg-3', component_property='value'),
        Input(component_id='color-by-feature-3', component_property='value')
    ]
)
def update_dim_reduc_3(tumor, algo, gene):
    return _build_dim_reduc(tumor, algo, gene, 'gene')


@app.callback(
    Output(component_id='dim-reduc-scatter-4', component_property='figure'),
    [
        Input(component_id='select-tumor-4', component_property='value'),
        Input(component_id='dim-reduc-alg-4', component_property='value'),
        Input(component_id='color-by-feature-4', component_property='value')
    ]
)
def update_dim_reduc_4(tumor, algo, gene):
    return _build_dim_reduc(tumor, algo, gene, 'gene')

def build_dim_reduc_selector(idd):
    return dcc.RadioItems(
        options=[
            {'label': 'UMAP', 'value': 'umap'},
            {'label': 'PHATE', 'value': 'phate'}
        ],
        value='umap',
        id=idd
    )

def build_num_dims_selector(idd):
    return dcc.RadioItems(
        options=[
            {'label': '2', 'value': 2},
            {'label': '3', 'value': 3}
        ],
        value=3,
        id=idd
    )

def build_feature_category_selector(idd):
    return dcc.RadioItems(
        options=[
            {'label': 'Gene', 'value': 'gene'},
            {'label': 'Cluster', 'value': 'cluster'}
        ],
        value='gene',
        id=idd
    )

def build_features_selector(idd, tumor, category):
    if category == 'gene':
        return dcc.Input(id=idd, value='OLIG1')
    elif category == 'cluster':
        return dcc.Dropdown(
            options=[{'label': 'Cluster', 'value': 'Cluster'}],
            value='Cluster',
            id=idd
        )

def _build_dim_reduc(tumor_id, algo, feat, category):
    if algo == 'umap':
        df_dim_reduc = load_data.load_tumor_umap(tumor_id, 3)
    elif algo == 'phate':
        df_dim_reduc = load_data.load_tumor_phate(tumor_id, 3)

    if category == 'gene':
        if feat in load_data.load_tumor_gene_names(tumor_id):
            df_color = load_data.load_tumor_gene(tumor_id, feat)
            col = 'color_by'
            color_range = [
                min(df_color[col]),
                max(df_color[col])
            ]
            df = df_dim_reduc.join(df_color)
            markers=dict(
                size=5,
                color=df[col],
                colorscale='Viridis',
                opacity=0.0,
                cmin=color_range[0],
                cmax=color_range[1],
                colorbar=dict(
                    thickness=20
                )
            )
            fig = go.Figure(data=[go.Scatter3d(
                x=df[df_dim_reduc.columns[0]],
                y=df[df_dim_reduc.columns[1]],
                z=df[df_dim_reduc.columns[2]],
                mode='markers',
                marker=markers,
                showlegend=False
            )])
    elif category == 'cluster':
        df_color = load_data.load_tumor_clusters_for_cells(tumor_id)
        col = 'color_by'
        df = df_dim_reduc.join(df_color)
        fig = go.Figure()
        for clust_i, clust in enumerate(sorted(set(df[col]))):
            df_clust = df.loc[df[col] == clust]
            markers=dict(
                size=5,
                color=PALETTE[clust_i],
                opacity=1.0
            )
            fig.add_trace(
                go.Scatter3d(
                    x=df_clust[df_clust.columns[0]],
                    y=df_clust[df_clust.columns[1]],
                    z=df_clust[df_clust.columns[2]],
                    mode='markers',
                    marker=markers,
                    name="Cluster {}".format(clust)
                )
            )
    fig.update_layout(
        autosize=True,
        width=FIG_DIM+10,
        height=FIG_DIM
    )
    return fig

def _build_reduc_card(title, graph_id):
    return dbc.Card(
        [
            dbc.CardHeader(
                title,
                style={
                    "background-color":"#e3e3e3",
                    "font-weight":"bold",
                    "font-size":"Large",
                    "text-align": "center"
                }
            ),
            dbc.CardBody([
                dcc.Graph(
                    id=graph_id,
                    figure=_build_dim_reduc('PJ017', 'umap', 'OLIG1', 'gene')
                )
            ], style={"height": "10vh", "width": "100%"})
        ],
        style={"height": SCATTER_HEIGHT, "width": "100%"}
        #style={"height": "100%", "width": "50%"}
    )


def _build_control_panel():
    def _build_control_panel_for_plot(plot_num):
        return [
            html.H5(
                "Plot {}".format(plot_num),
                style={
                    'text-align': 'center',
                    "font-weight":"bold"
                }
            ),
            html.H6("Select a tumor to visualize:"),
            common.build_tumor_dropdown('select-tumor-{}'.format(plot_num)),
            html.H6("Select dimension reduction:"),
            build_dim_reduc_selector('dim-reduc-alg-{}'.format(plot_num)),
            html.H6("Select dimensions:"),
            build_num_dims_selector('num-dims-{}'.format(plot_num)),
            html.H6("Select dimensions:"),
            build_feature_category_selector('select-feature-category-{}'.format(plot_num)),
            html.H6("Color by feature (enter a gene): "),
            #dcc.Input(id='color-by-gene-{}'.format(plot_num), value='OLIG1'),
            html.Div([
                build_features_selector('color-by-feature-{}'.format(plot_num), 'PJ016', 'gene') 
            ], id='color-by-feature-container-{}'.format(plot_num))
        ]
    control_panel = []
    for i in [1,2,3,4]:
        control_panel += _build_control_panel_for_plot(i)
        if i < 4:
            control_panel.append(html.Hr())
    return control_panel


LAYOUT = dcc.Tab(
    label='Dimension Reduction',
    children=[
        dbc.Container(fluid=True, children=[
            dbc.Row(html.Hr(), style={'height': '1%'}),
            dbc.Row(children=[
                dbc.Col(width=100, style={'width': '1%'}),
                dbc.Col(
                    [
                        dbc.Card(children=[
                            dbc.CardHeader(
                                "Control Panel",
                                style={
                                    "background-color":"#e3e3e3",
                                    "font-weight":"bold",
                                    "font-size":"Large",
                                    "text-align": "center"
                                }
                            ),
                            dbc.CardBody(
                                _build_control_panel()
                            )
                        ])
                    ],
                    width='100',
                    style={'width': '15%'}
                ),
                dbc.Col(width=100, style={'width': '1%'}),
                dbc.Row(
                    [
                        dbc.Col(
                            [
                                 _build_reduc_card('Plot 1', 'dim-reduc-scatter-1'),
                                 dbc.Row(html.Hr(), style={'height': '1%'}),
                                 _build_reduc_card('Plot 3', 'dim-reduc-scatter-3')
                            ],
                            width={"size": "auto", "order": 1}
                        ),
                        dbc.Col(
                            [
                                 _build_reduc_card('Plot 2', 'dim-reduc-scatter-2'),
                                 dbc.Row(html.Hr(), style={'height': '1%'}),
                                 _build_reduc_card('Plot 4', 'dim-reduc-scatter-4')
                            ],
                            width={"size": "100px", "order": 2}
                        )
                    ], 
                    style={"width": "80%"}
                )
            ])
        ])
    ]
)
