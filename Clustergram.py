import pandas as pd

import dash
import dash_bio as dashbio
import dash_html_components as html
import dash_core_components as dcc

correlation_matrix = pd.read_csv("Data/GTEx/CorrelationMatrix/corr_matrix.tsv", 
    header=0, index_col=0, sep="\t")

gtex = pd.read_csv("Data/GTEx/Metadata/GTEx.tsv", 
    header=0, index_col=0, sep="\t")

rows = list(correlation_matrix.index)
columns = list(correlation_matrix.columns.values)

external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div([
    "Rows to display",
    dcc.Dropdown(
        id='clustergram-input',
        value=rows[:],
        multi=True,
        options=[{'label': row, 'value': row} \
                for row in list(correlation_matrix.index)]),

    html.Div(id='my-clustergram')
])

@app.callback(
    dash.dependencies.Output('my-clustergram', 'children'),
    [dash.dependencies.Input('clustergram-input', 'value')]
)
def update_clustergram(rows):
    if len(rows) < 2:
        return "Please select at least two rows to display."

    return dcc.Graph(figure=dashbio.Clustergram(
        data=correlation_matrix.loc[rows].values,
        column_labels=columns,
        row_labels=rows,
        color_threshold={'row': 250, 'col': 700},
        hidden_labels=['col', 'row'],
        height=800, width=1000,
        optimal_leaf_order=True,
        color_map=[[0.0, '#636EFA'], [0.25, '#AB63FA'],
                [0.5, '#FFFFFF'], [0.75, '#E763FA'], [1.0, '#EF553B']]
        ))
        
if __name__ == '__main__':
    app.run_server(debug=True)