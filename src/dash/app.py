import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import requests

# Fontawesome (icons)
FA = "https://use.fontawesome.com/releases/v5.8.1/css/all.css"

# Logo
LOGO = "https://picsum.photos/30/30"

# App
app = dash.Dash(__name__,
                external_stylesheets=[dbc.themes.BOOTSTRAP, FA],
                prevent_initial_callbacks=True,
                suppress_callback_exceptions=True)

### Sidebar

# Menu 1
menu_1 = [
    html.Li(
        dbc.Row(
            [
                dbc.Col(    html.I(className="fas fa-table mr-3"), width=1, className="icon-padding" ),
                dbc.Col(    dbc.NavLink("Page 1", href="/page-1/"), className="navlink-padding"  )
            ],
            className="my-1",
        ),
        id="menu-1",
    )
]

# Menu 2
menu_2 = [
    html.Li(
        dbc.Row(
            [
                dbc.Col(    html.I(className="fas fa-chart-line mr-3"), width=1, className="icon-padding"  ),
                dbc.Col(    dbc.NavLink("Page 2", href="/page-2/"), className="navlink-padding"  ),
            ],
            className="my-1",
        ),
        id="menu-2",
    )
]

# Sidebar
sidebar = html.Div(
    [
        dbc.Nav(menu_1 + menu_2, vertical=True),
    ],
    className="sidebar",
    id="sidebar",
)

### Nav Bar

# Nav Bar Row
navbarRow = dbc.Row(
    [
    ],
    no_gutters=True,
    className="ml-auto flex-nowrap mt-3 mt-md-0",
    align="center",
)

# Nav Bar
navbar = dbc.Navbar(
    [
        html.A(
            dbc.Row(
                [
                    dbc.Col(html.Img(src=LOGO, height="30px")),
                    dbc.Col(dbc.NavbarBrand("Dash Template", className="ml-2")),
                ],
                align="center",
                no_gutters=True,
            ),
            href="https://github.com/ckrusemd",
        ),
        dbc.NavbarToggler(id="navbar-toggler"),
        dbc.Collapse(navbarRow, id="navbar-collapse", navbar=True),
    ],
    color="dark",
    dark=True,
)

# Content
content = html.Div(id="page-content", className="content")

# Layout
app.layout = html.Div([dcc.Location(id="url"), navbar,sidebar, content])

### Page Layouts

# Page 1
page_1 = html.Div([
    html.H1("Page 1"),
    html.Button('Test API', id='test-api', n_clicks=0),
    html.Div(id='container-api-result',
             children='')
])

# Page 2
page_2 = html.Div([
    html.H1("Page 2")
])

# URL
@app.callback(Output("page-content", "children"), [Input("url", "pathname")])
def render_page_content(pathname):
    if pathname in ["/", "/page-1/"]:
        return page_1
    elif pathname == "/page-2/":
        return page_2
    return dbc.Jumbotron(
        [
            html.H1("404: Not found", className="text-danger"),
            html.Hr(),
            html.P(f"The pathname {pathname} was not recognised..."),
        ]
    )

# Test the API
@app.callback(
    dash.dependencies.Output('container-api-result', 'children'),
    [dash.dependencies.Input('test-api', 'n_clicks')])
def update_output(n_clicks):
    greeting = requests.get("http://127.0.0.1:80/test/HelloWorld").json()
    return html.H2(greeting)

# Run the server
if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')