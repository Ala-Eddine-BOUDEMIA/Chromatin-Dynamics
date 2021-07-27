import numpy as np
import pandas as pd
import networkx as nx
from pyvis import network as net

import Config

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

def draw_graph3(
	networkx_graph, notebook = True,
	output_filename = 'Unweighted_Networks/TCGA/graph.html',
	show_buttons = True,
	only_physics_buttons = False):

    pyvis_graph = net.Network(notebook = notebook)
    pyvis_graph.width = '2000px'

    for node, node_attrs in networkx_graph.nodes(data = True):
        pyvis_graph.add_node(node, title = str(node), **node_attrs)

    for source,target,edge_attrs in networkx_graph.edges(data = True):
        if not 'value' in edge_attrs and not 'width' in edge_attrs and 'weight' in edge_attrs:
            edge_attrs['value'] = edge_attrs['weight']
        pyvis_graph.add_edge(source,target,**edge_attrs)

    if show_buttons:
        if only_physics_buttons:
            pyvis_graph.show_buttons(filter_ = ['physics'])
        else:
            pyvis_graph.show_buttons()

    return pyvis_graph.show(output_filename)

def adjacencyMatrix(
	cv_list, g_corr):

	correlation_matrix = pd.read_csv(g_corr, 
		header = 0, index_col = 0, sep = '\t')

	cv_list = pd.read_csv(cv_list, 
		header = 0, index_col = 0, sep = ';')
	
	for i in range(len(correlation_matrix.index)):
		for j in range(len(correlation_matrix.columns)):
			if i == j:
				correlation_matrix.iloc[i, j] = 0

	correlation_matrix = correlation_matrix.join(cv_list["GeneName"])
	labels = correlation_matrix["GeneName"]
	
	adjacency_matrix = pd.DataFrame(np.where(
		correlation_matrix.iloc[:, correlation_matrix.columns != "GeneName"] > 0.7, 1, 0), 
		index = labels,
		columns = labels)
	
	G = nx.from_pandas_adjacency(adjacency_matrix)

	draw_graph3(G)

if __name__ == '__main__':

	adjacencyMatrix(
		cv_list = Config.args.list,
		g_corr = Config.args.corrFullG)