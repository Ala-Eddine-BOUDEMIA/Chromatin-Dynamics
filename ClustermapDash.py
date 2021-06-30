import pandas as pd

import dash
import dash_bio as dashbio
import dash_html_components as html
import dash_core_components as dcc

import Config 

counts = pd.read_csv(Config.args.cv, 
    header = 0, index_col = 0, sep = "\t")

metadata = pd.read_csv(Config.args.meta, 
    header = 0, index_col = 0, sep = "\t")

cv_list = pd.read_csv(Config.args.list,
    header = 0, index_col = 0, sep = ";")

counts = counts.join(cv_list["GeneName"])

rows = list(counts["GeneName"])
columns = list(counts.columns.values)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    "Rows to display",
    dcc.Dropdown(
        id = 'clustergram-input',
        value = rows[:],
        multi = True,
        options = [{'label': row, 'value': row} \
                for row in list(counts["GeneName"])]),

    html.Div(id = 'my-clustergram')
])

@app.callback(
    dash.dependencies.Output('my-clustergram', 'children'),
    [dash.dependencies.Input('clustergram-input', 'value')]
)
def update_clustergram(rows):
    if len(rows) < 2:
        return "Please select at least two rows to display."

    return dcc.Graph(
        figure = dashbio.Clustergram(
        data = counts.iloc[:, 8000:9000].values,
        column_labels = columns,
        row_labels = rows,
        color_threshold = {'row': 250, 'col': 700},
        hidden_labels = ['col'],
        height = 1200, width = 1400,
        optimal_leaf_order = True,
        color_map = [[0.0, '#636EFA'], [0.25, '#AB63FA'],
                [0.5, '#FFFFFF'], [0.75, '#E763FA'], [1.0, '#EF553B']]
        ))
        
if __name__ == '__main__':
    app.run_server(debug=True)