import pandas as pd
import plotly.express as px

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

gtex_counts = pd.read_csv("Test/Data/GTEx/read_counts.tsv", 
					header=0, index_col=0, sep="\t")

gtex = pd.read_csv("Test/Data/GTEx/GTEx.tsv", 
					header=0, index_col=0, sep="\t")

mean = gtex_counts.mean(axis=1)
std = gtex_counts.std(axis=1)

df = pd.DataFrame(data=[mean, std], index=['mean','std']).T

print(df.columns)

fig = px.scatter(df, x="mean", y="std")
fig.show()

mean = gtex_counts.mean(axis=0)
std = gtex_counts.std(axis=0)

df = pd.DataFrame(data=[mean, std], index=['mean','std']).T

fig = px.scatter(df, x="mean", y="std")
fig.show()
