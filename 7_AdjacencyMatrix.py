import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import plotly.graph_objects as go

import Config

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', 10)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', 1)

def adjacencyMatrix(
	meta, cv_list, g_corr_cv):

	correlation_matrix = pd.read_csv(g_corr_cv, 
		header = 0, index_col = 0, sep = '\t')

	cv_list = pd.read_csv(cv_list, 
		header = 0, index_col = 0, sep = '\t')

	correlation_matrix = correlation_matrix.join(cv_list["GeneName"])
	labels = correlation_matrix["GeneName"]

	adjacency_matrix = pd.DataFrame(np.where(
		correlation_matrix.iloc[:, correlation_matrix.columns != "GeneName"] > 0.7, 1, 0), 
		index = labels,
		columns = labels)

	G = nx.from_pandas_adjacency(adjacency_matrix)
	G.name = "Graph from pandas adjacency matrix"
	nx.draw_kamada_kawai(G, 
		with_labels = True, 
		edge_cmap = "inferno")
	
	plt.show()


if __name__ == '__main__':

	adjacencyMatrix(
		meta = Config.args.meta,
		cv_list = Config.args.list,
		g_corr_cv = Config.args.corrCVg)