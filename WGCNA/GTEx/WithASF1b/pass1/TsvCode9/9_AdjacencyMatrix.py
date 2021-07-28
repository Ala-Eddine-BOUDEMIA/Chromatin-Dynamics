import numpy as np
import pandas as pd
import networkx as nx
from pyvis import network as net

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

def draw_graph3(
	networkx_graph, notebook = True,
	output_filename = 'graph.html',
	show_buttons = True,
	only_physics_buttons = False):
	
	membership = pd.read_csv('membership.txt', header = 0, index_col = 0, sep = "\t")
	pyvis_graph = net.Network(notebook = notebook)
	pyvis_graph.width, pyvis_graph.height = '2000px', '1000px'
	
	for node, node_attrs in networkx_graph.nodes(data = True):
		pyvis_graph.add_node(node, title = str(node), group = membership.loc[node]["Module"], **node_attrs)
		
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

def adjacencyMatrix():

	adjacency_matrix = pd.read_csv('pass1.tsv', header = 0, index_col = 0, sep = "\t")
	G = nx.from_pandas_adjacency(adjacency_matrix)

	draw_graph3(G,)

if __name__ == '__main__':

	adjacencyMatrix()