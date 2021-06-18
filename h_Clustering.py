import numpy as np
import pandas as pd
import seaborn as  sns

from matplotlib.patches import Patch
from scipy.spatial.distance import pdist, squareform

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

correlation_matrix = correlation_matrix.join(gtex["smts"])
correlation_matrix = correlation_matrix.dropna()
tissues = correlation_matrix.pop("smts")

palette1 = sns.hls_palette(10)
palette2 = sns.color_palette("bwr",10)
palette3 = sns.color_palette("inferno",10)
palette = palette1 + palette2 + palette3
lut = dict(zip(set(tissues.unique()), palette))
colors = tissues.map(lut)

g = sns.clustermap(correlation_matrix, 
                vmin=0, vmax=1, 
                cmap="icefire", standard_scale=1,
                row_colors=colors, col_colors=colors, 
                xticklabels=False, yticklabels=False,
                method="single")

handles = [Patch(facecolor=lut[name]) for name in lut]
g.ax_row_dendrogram.legend(handles, lut, title='Tissues',
        bbox_to_anchor=(0, 1), loc='best')

g.savefig("Images/GTEx/Clustermap/Clustergram_full.png", dpi = 300)