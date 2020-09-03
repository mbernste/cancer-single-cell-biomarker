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
import nav


def charts():
    layout = html.Div(children=[
        nav.LAYOUT,
        html.H1(children='CHARTS: CHARacterizing Tumor Subpopulations'),

        html.Div(children='''
            Dash: A web application framework for Python.
        '''),

        dbc.Container(fluid=True, children=[
            dcc.Tabs([
                dim_reduc.dim_reduc(),
                #clust_compare.LAYOUT,
            ])
        ])
    ])
    return layout
