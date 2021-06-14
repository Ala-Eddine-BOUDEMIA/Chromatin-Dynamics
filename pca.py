import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

gtex_counts = pd.read_csv("Data/GTEx/Normalized/GTExNormalized.tsv", 
					header=0, index_col=0, sep="\t")

gtex = pd.read_csv("Data/GTEx/GTEx.tsv", 
					header=0, index_col=0, sep="\t")

tissues = gtex["smts"]

# code for testing
"""
tissues_set = ["Bladder", "Bone Marrow", "Cervix Uteri", "Fallopian Tube"]

for index in tissues.index:
	if tissues.loc[index] not in tissues_set:
		tissues.drop(index, axis=0, inplace=True)
"""

# top 1000 genes
"""
gtex_counts["total"] = gtex_counts.sum(axis=0)
gtex_counts.sort_values("total", ascending=False, inplace = True)
gtex_counts = pd.DataFrame(gtex_counts.iloc[1:1000]) 
gtex_counts.to_csv("gtex_counts_top1000.tsv", sep="\t")
"""

pca = PCA().fit(gtex_counts.dropna().T)
Y = pca.fit_transform(gtex_counts.dropna().T)
ratio = pca.explained_variance_ratio_ * 100

d = pd.DataFrame(index=gtex_counts.T.index)
d["PC1"] = Y[:, 0]
d["PC2"] = Y[:, 1]
d = d.join(tissues)

d.to_csv("pca_top1000.tsv", sep="\t")

print("PC1: ", round(ratio[0],2))
print("PC2: ", round(ratio[1],2))

d = pd.read_csv("pca_top1000.tsv", 
				header=0, index_col=0, sep="\t")

fig = px.scatter(
    d.dropna(), x="PC1", y="PC2",
    color="smts",
    title="GTEx PCA")
fig.show()