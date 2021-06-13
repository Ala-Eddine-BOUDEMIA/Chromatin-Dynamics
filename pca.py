import pandas as pd
import plotly.express as px

from sklearn.decomposition import PCA

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

gtex_counts = pd.read_csv("Test/Data/GTEx/read_counts.tsv", 
					header=0, index_col=0, sep="\t")

gtex = pd.read_csv("Test/Data/GTEx/GTEx.tsv", 
					header=0, index_col=0, sep="\t")

tissues = gtex["smts"]
tissues_set = ["Bladder", "Bone Marrow", "Cervix Uteri", "Fallopian Tube"]

for index in tissues.index:
	if tissues.loc[index] not in tissues_set:
		tissues.drop(index, axis=0, inplace=True)

pca = PCA().fit(gtex_counts.dropna().T)
Y = pca.fit_transform(gtex_counts.dropna().T)
ratio = pca.explained_variance_ratio_ * 100

d = pd.DataFrame(index=gtex_counts.T.index)
d["PC1"] = Y[:, 0]
d["PC2"] = Y[:, 1]
d = d.join(tissues)

fig = px.scatter(
    d, x="PC1", y="PC2",
    color="smts",
    title="GTEx PCA")
fig.show()