import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.figure_factory as ff

from scipy.spatial.distance import pdist, squareform
import seaborn as  sns

# metadata
gtex = pd.read_csv("Data/GTEx/Metadata/GTEx.tsv", 
	header=0, index_col=0, sep="\t")

# complete dataset
"""
gtex_counts = pd.read_csv("Data/GTEx/Normalized/GTExNormalized.tsv", 
	header=0, index_col=0, sep="\t")
"""
# top 1000 genes
"""
gtex_counts = pd.read_csv("Data/GTEx/Top1000/GTExTop1000Genes.tsv", 
	header=0, index_col=0, sep="\t")

gtex_counts = pd.DataFrame(np.log2(gtex_counts + 1))

correlation_matrix = gtex_counts.corr(method="pearson")
correlation_matrix.to_csv('Data/GTEx/CorrelationMatrix/corr_matrix.tsv', sep="\t")
"""
# correlation matrix
correlation_matrix = pd.read_csv("Data/GTEx/CorrelationMatrix/corr_matrix.tsv",
	header=0, index_col=0, sep="\t")

correlation_matrix = correlation_matrix.join(gtex["smtsd"])
tissues = correlation_matrix.pop("smtsd")

lut = dict(zip(set(tissues.unique()), sns.hls_palette(len(set(tissues)))))
colors = tissues.map(lut)

g=sns.clustermap(correlation_matrix, vmin=-1, vmax=1, 
                 row_colors=colors, col_colors=colors, 
                 xticklabels=False, yticklabels=False,
                 method="complete")

g.savefig("Images/GTEx/Clustermap/Clustergram_full.png", dpi = 300)