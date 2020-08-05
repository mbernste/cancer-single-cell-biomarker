import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_bootstrap_components as dbc

def build_tumor_dropdown(html_id):
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
            {'label': 'PJ018 + PJ025', 'value': 'PJ018&PJ025'},
        ],
        value='PJ016',
        id=html_id
    )
