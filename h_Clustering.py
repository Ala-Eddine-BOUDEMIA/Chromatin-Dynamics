import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.figure_factory as ff

from scipy.spatial.distance import pdist, squareform
import seaborn as  sns

# complete dataset
"""
gtex_counts = pd.read_csv("Data/GTEx/Normalized/GTExNormalized.tsv", 
	header=0, index_col=0, sep="\t")
"""
# top 1000 genes
"""
gtex_counts = pd.read_csv("Data/GTEx/Top1000/GTExTop1000Genes.tsv", 
	header=0, index_col=0, sep="\t")

gtex_counts = pd.DataFrame(np.log2(gtex_counts + 1))

correlation_matrix = gtex_counts.corr(method="pearson")
correlation_matrix.to_csv('corr_matrix.tsv', sep="\t")
"""

correlation_matrix = pd.read_csv("1000vs1000.tsv", 
	header=0, index_col=0, sep="\t")

gtex = pd.read_csv("Data/GTEx/Metadata/GTEx.tsv", 
	header=0, index_col=0, sep="\t")

tissues = gtex["smtsd"]


print("correlation matrix loaded")

# Initialize figure by creating upper dendrogram
fig = ff.create_dendrogram(correlation_matrix, 
    orientation='bottom')

for i in range(len(fig['data'])):
    fig['data'][i]['yaxis'] = 'y2'

print("Upper dendrogram created")

# Create Side Dendrogram
dendro_side = ff.create_dendrogram(correlation_matrix, 
    orientation='right')

for i in range(len(dendro_side['data'])):
    dendro_side['data'][i]['xaxis'] = 'x2'

print("Side dendrogram created")

# Add Side Dendrogram Data to Figure
for data in dendro_side['data']:
    fig.add_trace(data)

# Create Heatmap
dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
dendro_leaves = list(map(int, dendro_leaves))
data_dist = pdist(correlation_matrix.to_numpy())
heat_data = squareform(data_dist)
heat_data = heat_data[dendro_leaves,:]
heat_data = heat_data[:,dendro_leaves]

heatmap = [go.Heatmap(
    x = dendro_leaves,
    y = dendro_leaves,
    z = heat_data,
    colorscale = 'inferno',
    hoverinfo = "skip")]

heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

print("Heatmap created")

# Add Heatmap Data to Figure
for data in heatmap:
    fig.add_trace(data)

# Edit Layout
fig.update_layout({'width':1000, 'height':1000,
    'showlegend':False, 'hovermode': 'closest'})

# Edit xaxis
fig.update_layout(xaxis={'domain': [.15, 1],
    'mirror': False, 'showgrid': False, 'showline': False,
    'zeroline': False, 'ticks':""})

# Edit xaxis2
fig.update_layout(xaxis2={'domain': [0, .15],
    'mirror': False, 'showgrid': False, 'showline': False,
    'zeroline': False, 'showticklabels': False, 'ticks':""})

# Edit yaxis
fig.update_layout(yaxis={'domain': [0, .85],
    'mirror': False, 'showgrid': False, 'showline': False,
    'zeroline': False, 'showticklabels': False, 'ticks': ""})

# Edit yaxis2
fig.update_layout(yaxis2={'domain':[.825, .975],
    'mirror': False, 'showgrid': False, 'showline': False,
    'zeroline': False, 'showticklabels': False, 'ticks':""})

# Plot
fig.show()
fig.write_html("Plotly_HTML_Files/GTEx/H_Clustering/hcluster.html")