import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import requests
from loguru import logger

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
                dbc.Col(    dbc.NavLink("Create Project", href="/create_project/"), className="navlink-padding"  )
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
                dbc.Col(    dbc.NavLink("Current Project", href="/current_project/"), className="navlink-padding"  ),
            ],
            className="my-1",
        ),
        id="menu-2",
    )
]

# Menu 3
menu_3 = [
    html.Li(
        dbc.Row(
            [
                dbc.Col(    html.I(className="fas fa-chart-line mr-3"), width=1, className="icon-padding"  ),
                dbc.Col(    dbc.NavLink("Review Abstracts", href="/review_abstracts/"), className="navlink-padding"  ),
            ],
            className="my-1",
        ),
        id="menu-3",
    )
]

# Sidebar
sidebar = html.Div(
    [
        dbc.Nav(menu_1 + menu_2 + menu_3, vertical=True),
    ],
    className="sidebar",
    id="sidebar",
)

### Nav Bar

# Nav Bar Row
navbarRow = dbc.Row(
    [
        dbc.Col( dbc.NavbarBrand("Current Project",className="navbar-offset") ),
        dbc.Col( dcc.Dropdown(
            id='current-project',
            options=[
                {'label': 'Kruse Eiken Vestergaard', 'value': 'Kruse Eiken Vestergaard'}
            ],
            value='Kruse Eiken Vestergaard'
        ),width=6),
        dbc.Col( dbc.NavbarBrand("Current User",className="navbar-offset") ),
        dbc.Col(dcc.Dropdown(
            id='demo-dropdown',
            options=[
                {'label': 'Christian Kruse', 'value': 'Christian Kruse'}
            ],
            value='Christian Kruse'
        ),width=6)
    ],
    no_gutters=True,
    className="ml-4 flex-nowrap",
    align="center",
)

# Nav Bar
navbar = dbc.Navbar(
    [
        dbc.Col(
            html.A(
                dbc.Row(
                    [
                        dbc.Col(html.Img(src=LOGO, height="30px")),
                        dbc.Col(dbc.NavbarBrand("Meta-Analysis Tool", className="ml-2 flex-nowrap")),
                        
                    ],
                    align="center",
                    no_gutters=True,
                ),
                href="https://github.com/ckrusemd",
            ),width=2),
        dbc.Col(width=4),
        dbc.Col(
            dbc.Row([
                dbc.Collapse(navbarRow, id="navbar-collapse", navbar=True)
            ]),width=6),
    ],
    color="dark",
    dark=True,
)

# Content
content = html.Div(id="page-content", className="content")

# Layout
app.layout = html.Div([dcc.Location(id="url"), navbar,sidebar, content])

### Page Layouts

# Create Project
page_1 = html.Div([
    html.H1("Create Project"),
    html.H5("Project Name:"),
    dcc.Input( id="project_Name" ,placeholder="hello",value="Kruse Eiken Vestergaard"),
    html.H5("Project Description:"),
    dcc.Textarea( id="project_Description", style={'width': '100%', 'height': 300} ,value="My publications"),
    html.H5("Project Query:"),
    dcc.Input( id="project_Query" ,value="Kruse Eiken Vestergaard"),
    html.H5("Project Author:"),
    dcc.Input( id="project_Author" ,value=1),
    html.P(),
    html.Button('Save Project', id='save-project-button', n_clicks=0),
    html.Div(id="create-project-content"),
    html.Div(id='container-api-result',
             children='')
])

# Current Project
page_2 = html.Div([
    html.H1("Current Project"),
    html.Button("Get Current Project", id='get-current-project', n_clicks=0),
    html.Div(id="current-project-content")
])

# Review Abstracts
page_3 = html.Div([
    html.H1("Review Abstracts"),
    html.Button("Get Next Abstract", id='get-next-abstract', n_clicks=0)
])

# Get Current Project
@app.callback(
    dash.dependencies.Output('current-project-content', 'children'),
    [dash.dependencies.Input('get-current-project', 'n_clicks')])
def get_current_project(n_clicks):
    projects = requests.get("http://api:8080/mongodb_projects/projects").json()
    logger.info( projects )
    result = []
    for project in projects:
        for key,value in project.items():
            result = result + [html.P(key + ": " + str( value ) )]
    return result

# Create New Project
@app.callback(
    dash.dependencies.Output('create-project-content', 'children'),
    [   dash.dependencies.Input('save-project-button', 'n_clicks'),
        dash.dependencies.State('project_Name', 'value'),
        dash.dependencies.State('project_Description', 'value'),
        dash.dependencies.State('project_Query', 'value'),
        dash.dependencies.State('project_Author', 'value')
    ])
def create_new_project(n_clicks,project_Name,project_Description,project_Query,project_Author):
    post_body = {"name": project_Name,"description": project_Description,"query": project_Query,"user_id": project_Author}
    result = requests.post("http://api:8080/mongodb_projects/create_project/",json=post_body).json()
    logger.info( result )
    return html.P("OK!")



    dcc.Input( id="project_Name" ,placeholder="hello"),
    html.H5("Project Description:"),
    dcc.Textarea( id="project_Description", style={'width': '100%', 'height': 300} ),
    html.H5("Project Query:"),
    dcc.Input( id="project_Query" ),
    html.H5("Project Author:"),
    dcc.Input( id="project_Author" ),
    html.P(),
    html.Button('Save Project', id='save-project-button', n_clicks=0),
    html.Div(id='container-api-result',
             children='')














# URL
@app.callback(Output("page-content", "children"), [Input("url", "pathname")])
def render_page_content(pathname):
    if pathname in ["/", "/create_project/"]:
        return page_1
    elif pathname == "/current_project/":
        return page_2
    elif pathname == "/review_abstracts/":
        return page_3


        
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
    greeting = requests.get("http://api:8080/test/HelloWorld").json()
    return html.H2(greeting)

# Run the server
if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0')