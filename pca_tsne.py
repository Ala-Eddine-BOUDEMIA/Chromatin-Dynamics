import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.manifold import TSNE 
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# complete dataset
gtex_counts = pd.read_csv("Data/GTEx/v_c/hc_hv.tsv", 
					header=0, index_col=0, sep="\t")

genes = gtex_counts.pop("Gene Official Symbol")

# top 1000 genes
"""
gtex_counts = pd.read_csv("Data/GTEx/Top1000/GTExTop1000Genes.tsv", 
					header=0, index_col=0, sep="\t")
"""
gtex_counts = pd.DataFrame(np.log2(gtex_counts + 1))

gtex = pd.read_csv("Data/GTEx/Metadata/GTEx.tsv", 
					header=0, index_col=0, sep="\t")


tissues = gtex["smtsd"]
lib_size = gtex_counts.sum(axis=0).to_frame(name="lib_size")

scaler = StandardScaler()
std_counts = scaler.fit_transform(gtex_counts.dropna().T)

pca = PCA(n_components=2)
P = pca.fit_transform(std_counts)
ratio = pca.explained_variance_ratio_ * 100

d = pd.DataFrame(index=gtex_counts.T.index)
d["PC1"] = P[:, 0]
d["PC2"] = P[:, 1]
d = d.join(lib_size)
d = d.join(tissues)
d.sort_values("smtsd", inplace=True)

print("PC1: ", round(ratio[0],2))
print("PC2: ", round(ratio[1],2))

"""
d.to_csv("Data/GTEx/PCA/pca_full.tsv", sep="\t")

d = pd.read_csv("Data/GTEx/PCA/pca.tsv", 
	header=0, index_col=0, sep="\t")
"""

fig = px.scatter(
    d.dropna(), 
    x="PC1", y="PC2",
    color="smtsd", 
    hover_data=[d.dropna().index, "lib_size"],
    #size="lib_size",
    title="GTEx PCA")
fig.show()
fig.write_html("Plotly_HTML_files/GTEx/PCA/hc_hv.html")

# Leading genes
loading_scores = pd.Series(pca.components_[0], index=gtex_counts.index)
sorted_loading_scores = loading_scores.abs().sort_values(ascending=False)
top_10_genes = sorted_loading_scores[0:10].index.values
print(loading_scores[top_10_genes])

# T-sne
tsne = TSNE(n_components=2)
T = tsne.fit_transform(std_counts)
df = pd.DataFrame(index=gtex_counts.T.index)
df["T1"] = T[:, 0]
df["T2"] = T[:, 1]
df = df.join(lib_size)
df = df.join(tissues)
df.sort_values("smtsd", inplace=True)

"""
df.to_csv("Data/GTEx/T-Sne/t_sne_full.tsv", sep="\t")

d = pd.read_csv("Data/GTEx/T-Sne/t_sne_full.tsv", 
	header=0, index_col=0, sep="\t")
"""

fig = px.scatter(
    df.dropna(), 
    x="T1", y="T2",
    color="smtsd", 
    hover_data=[df.dropna().index, "lib_size"],
    #size="lib_size",
    title="GTEx t-sne")
fig.show()
fig.write_html("Plotly_HTML_files/GTEx/T-Sne/hc_hv.html")