import pandas as pd
import numpy as np
import plotly.express as px

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

gtex_counts = pd.read_csv("Data/GTEx/Normalized/GTExNormalized.tsv", 
					header=0, index_col=0, sep="\t")

gtex = pd.read_csv("Test/Data/GTEx/GTEx.tsv", 
					header=0, index_col=0, sep="\t")

# Mean-Variance Gene
# f(std) = mean
mean = gtex_counts.mean(axis=1)
std = gtex_counts.std(axis=1)
df = pd.DataFrame(data=[mean, std], index=['mean','std']).T

fig = px.scatter(df, x="mean", y="std")
fig.show()

# f(std) = log2(mean+1)
mean_log10 = np.log10(mean+1) 
df_mean_log2 = pd.DataFrame(data=[mean_log10, std], index=['log10(mean+1)','std']).T

fig = px.scatter(df_mean_log2, x="log10(mean+1)", y="std")
fig.show()

# f(std) = mean(log2(gtex_counts+1))
mean_log10_counts = pd.DataFrame(np.log10(gtex_counts + 1)).mean(axis=1) 
df_mean_log10_counts = pd.DataFrame(data=[mean_log10_counts, std], index=['mean(log10(gtex_counts+1))','std']).T

fig = px.scatter(df_mean_log10_counts, x="mean(log10(gtex_counts+1))", y="std")
fig.show()

# Mean-Variance Sample
# f(std) = mean
# f(std) = mean
mean = gtex_counts.mean(axis=0)
std = gtex_counts.std(axis=0)
df = pd.DataFrame(data=[mean, std], index=['mean','std']).T

fig = px.scatter(df, x="mean", y="std")
fig.show()

# f(std) = log2(mean+1)
mean_log10 = np.log10(mean+1) 
df_mean_log2 = pd.DataFrame(data=[mean_log10, std], index=['log10(mean+1)','std']).T

fig = px.scatter(df_mean_log2, x="log10(mean+1)", y="std")
fig.show()

# f(std) = mean(log2(gtex_counts+1))
mean_log10_counts = pd.DataFrame(np.log10(gtex_counts + 1)).mean(axis=0) 
df_mean_log10_counts = pd.DataFrame(data=[mean_log10_counts, std], index=['mean(log10(gtex_counts+1))','std']).T

fig = px.scatter(df_mean_log10_counts, x="mean(log10(gtex_counts+1))", y="std")
fig.show()