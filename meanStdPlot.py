import numpy as np
import pandas as pd
import plotly.express as px

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# Complete dataset
gtex_counts = pd.read_csv("Data/GTEx/Normalized/GTExNormalized.tsv", 
						header=0, index_col=0, sep="\t")

# top 1000 genes
"""
gtex_counts = pd.read_csv("Data/GTEx/Top1000/GTExTop1000Genes.tsv", 
					header=0, index_col=0, sep="\t")
"""

# COUNTS
mean = np.log2(gtex_counts + 0.5).mean(axis=1)
std = np.sqrt(gtex_counts).std(axis=1)
df = pd.DataFrame(data=[mean, std], 
	index=['log2(mean+0.5)','sqrt(std)']).T

fig = px.scatter(df, 
	x="log2(mean+0.5)", 
	y="sqrt(std)",
	hover_data=[df.index],
	title="Mean Variance using counts")
fig.show()

# CPM
cpm = gtex_counts.div(gtex_counts.sum(axis=0).div(1e6))
mean_cpm = np.log2(cpm + 0.5).mean(axis=1)
std_cpm = np.sqrt(cpm).std(axis=1)
df_cpm = pd.DataFrame(data=[mean_cpm, std_cpm], 
	index=['log2(mean_cpm+0.5)','sqrt(std_cpm)']).T

fig = px.scatter(df_cpm, 
	x="log2(mean_cpm+0.5)", 
	y="sqrt(std_cpm)",
	hover_data=[df.index],
	title="Mean Variance using CPMs")
fig.show()

# Tissue
gtex = pd.read_csv("Data/GTEx/Metadata/GTEx.tsv", 
					header=0, index_col=0, sep="\t")

mean = np.log2(gtex_counts + 0.5).mean(axis=0)
std = np.sqrt(gtex_counts).std(axis=0)
df_tissue = pd.DataFrame(data=[mean, std, gtex['smtsd'].T], 
	index=['log2(mean+0.5)','sqrt(std)', 'smtsd']).T
df_tissue.sort_values("smtsd", inplace=True)

fig = px.scatter(df_tissue.dropna(),
	x="log2(mean+0.5)", 
	y="sqrt(std)",
	color="smtsd",
	hover_data=[df_tissue.dropna().index],
	title="Mean Variance using samples and colored by tissue type")
fig.show()