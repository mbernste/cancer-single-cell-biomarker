import dash_bootstrap_components as dbc
import dash_html_components as html


CHARTS_LOGO = "https://github.com/mbernste/cancer-single-cell-biomarker/raw/master/img/charts_logo.png"

LAYOUT = dbc.NavbarSimple(
    children=[
      dbc.NavItem(dbc.NavLink("About", href="/info")),
      dbc.DropdownMenu(
         nav=True,
         in_navbar=True,
         label="Menu",
         children=[
            dbc.DropdownMenuItem("Entry 1"),
            dbc.DropdownMenuItem("Entry 2"),
            dbc.DropdownMenuItem(divider=True),
            dbc.DropdownMenuItem("Entry 3"),
                  ],
              ),
            ],
    brand="CHARTS: Characterizing Tumor Subpopulations",
    brand_style={"font-size": "200%"},
    brand_href="/charts",
    sticky="top",
)

LAYOUT1 = dbc.NavbarSimple(
    children=[
        html.A(
            # Use row and col to control vertical alignment of logo / brand
            dbc.Row(
                [
                    dbc.Col(html.Img(src=CHARTS_LOGO, height="30px")),
                    dbc.Col(dbc.NavbarBrand("CHARTS\nCharacterizing tumor subpopulations", className="ml-2")),
                ],
                align="left",
                no_gutters=True,
            )
        ),
        dbc.NavItem(dbc.NavLink("Page 1", href="#")),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("More pages", header=True, style={"textcolor": "black"}),
                dbc.DropdownMenuItem("Page 2", href="#"),
                dbc.DropdownMenuItem("Page 3", href="#"),
            ],
            nav=True,
            in_navbar=True,
            label="More",
        ),
    ],
    #brand="CHARTS: Characterizing Tumor Subpopulations",
    #brand_href="#",
    color="dark",
    dark=True,
)
