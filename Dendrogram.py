import sys
import numpy as np
import pandas as pd
import seaborn as sns
import plotly.figure_factory as ff
import scipy.cluster.hierarchy as sch

import Config 

sys.setrecursionlimit(1000000)

chaperones_variants = pd.read_csv(Config.args.cv, 
    header = 0, index_col = 0, sep = "\t")

metadata = pd.read_csv(Config.args.meta, 
    header = 0, index_col = 0, sep = "\t")

genes = pd.read_csv(Config.args.list, 
    header = 0, index_col = 0, sep = ";")

chaperones_variants_log2 = pd.DataFrame(np.log2(chaperones_variants + 1))

chaperones_variants_log2 = chaperones_variants_log2.join(genes["GeneName"])
chaperones_variants_log2.set_index("GeneName", inplace = True)
chaperones_variants_log2 = chaperones_variants_log2.append(metadata["smtsd"])
chaperones_variants_log2.columns = chaperones_variants_log2.columns + " " + chaperones_variants_log2.loc['smtsd']
chaperones_variants_log2.drop("smtsd", axis = 0, inplace = True)

names = chaperones_variants_log2.columns
fig = ff.create_dendrogram(chaperones_variants_log2.T,
    orientation = 'left', labels = names,
    linkagefun = lambda chaperones_variants_log2: sch.linkage(chaperones_variants_log2.T, 
        method = "average", metric = "correlation"),)
fig.update_layout(width = 800, height = 1800)
fig.show()
