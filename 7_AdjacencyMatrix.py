import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

import Config

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

def adjacencyMatrix(
	cv_list, g_corr_cv):

	correlation_matrix = pd.read_csv(g_corr_cv, 
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
	print(nx.info(G))
	G.name = "Graph from pandas adjacency matrix"
	#pos=nx.graphviz_layout(G,prog='dot')
	nx.draw_random(G, 
		with_labels = True, 
		edge_cmap = "inferno")
	
	plt.show()

if __name__ == '__main__':

	adjacencyMatrix(
		cv_list = Config.args.list,
		g_corr_cv = Config.args.corrCVg)