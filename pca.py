import numpy as np
import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

# code for testing
"""
tissues_set = ["Bladder", "Bone Marrow", "Cervix Uteri", "Fallopian Tube"]

for index in tissues.index:
	if tissues.loc[index] not in tissues_set:
		tissues.drop(index, axis=0, inplace=True)
"""
"""
gtex_counts = pd.read_csv("Data/GTEx/Normalized/GTExNormalized.tsv", 
					header=0, index_col=0, sep="\t")
"""

# top 1000 genes
gtex_counts = pd.read_csv("Data/GTEx/Top1000/GTExTop1000Genes.tsv", 
					header=0, index_col=0, sep="\t")

gtex_counts = pd.DataFrame(np.log2(gtex_counts + 1))

gtex = pd.read_csv("Data/GTEx/Metadata/GTEx.tsv", 
					header=0, index_col=0, sep="\t")

tissues = gtex["smtsd"]
lib_size = gtex_counts.sum(axis=0).to_frame(name="lib_size")

scaler = StandardScaler()
std_counts = scaler.fit_transform(gtex_counts.dropna().T)

pca = PCA(n_components=3)
P = pca.fit_transform(std_counts)
ratio = pca.explained_variance_ratio_ * 100

d = pd.DataFrame(index=gtex_counts.T.index)
d["PC1"] = P[:, 0]
d["PC2"] = P[:, 1]
d = d.join(tissues)
d = d.join(lib_size)
d.sort_values("smtsd", inplace=True)

d.to_csv("Data/GTEx/PCA/pca.tsv", sep="\t")

print("PC1: ", round(ratio[0],2))
print("PC2: ", round(ratio[1],2))
print("PC3: ", round(ratio[2],2))
"""
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